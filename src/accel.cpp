/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <chrono>
#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

    void Accel::addMesh(Mesh *mesh) {
        if (m_mesh)
            throw NoriException("Accel: only a single mesh is supported!");
        m_mesh = mesh;
        m_bbox = m_mesh->getBoundingBox();
    }

    void Accel::build() {

        if (!m_mesh) return;

        mode = AccelMode::SAH;

        auto start = std::chrono::high_resolution_clock::now();

        m_tree.clear();

        auto root = AccelNode(m_mesh->getBoundingBox(), m_mesh->getTriangleCount());

        for (int i = 0; i < root.indices.size(); i++) {
            root.indices[i] = i;
        }

        //添加根节点
        m_tree.emplace_back(root);

        //辅助队列，存储八叉树节点的索引
        std::queue<uint32_t> q;
        q.push(0);

        //构建二叉树
        while (!q.empty()) {
            //层次遍历
            uint32_t n = q.size();
            for (uint32_t k = 0; k < n; k++) {

                //出队
                auto ptr = q.front();
                q.pop();

                if (mode == AccelMode::SAH) {
                    buildSAH(q, ptr);
                } else if (mode == AccelMode::OctTree) {
                    buildOctTree(q, ptr);
                }
            }
            m_depth_tree++;
        }

        auto over = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(over - start).count();

        std::cout << "[time to build the tree]: " << time << "ms" << std::endl;
        std::cout << "[max depth of the tree]: " << m_depth_tree << std::endl;
        std::cout << "[node count of the tree]: " << m_count_node << std::endl;
        std::cout << "[leaf count of the tree]: " << m_count_leaf << std::endl;
    }

    void Accel::buildOctTree(std::queue<uint32_t> &q, uint32_t nodeIndex) {
        if (m_tree[nodeIndex].indices.size() < COUNT_MIN)
            return;
        if (m_depth_tree > DEPTH_OCT_MAX)
            return;

        m_tree[nodeIndex].child = m_tree.size();
        Vector3f center = m_tree[nodeIndex].bbox.getCenter();
        for (uint32_t i = 0; i < 8; i++) {
            Vector3f corner = m_tree[nodeIndex].bbox.getCorner(i);
            Vector3f minPoint, maxPoint;
            for (uint32_t j = 0; j < 3; j++) {
                minPoint[j] = std::min(center[j], corner[j]);
                maxPoint[j] = std::max(center[j], corner[j]);
            }

            BoundingBox3f bbox_sub(minPoint, maxPoint);
            AccelNode node_sub(bbox_sub);
            for (auto face: m_tree[nodeIndex].indices)
                if (bbox_sub.overlaps(m_mesh->getBoundingBox(face)))
                    node_sub.indices.emplace_back(face);

            q.push(m_tree.size());
            m_tree.emplace_back(node_sub);
        }
        m_count_node += 8;
        m_count_leaf += 7;
    }

    void Accel::buildSAH(std::queue<uint32_t> &q, uint32_t nodeIndex) {
        if (m_tree[nodeIndex].indices.size() < COUNT_MIN)
            return;
        if (m_depth_tree > DEPTH_BVH_MAX)
            return;

        m_tree[nodeIndex].child = m_tree.size();
        uint32_t dimension = m_tree[nodeIndex].bbox.getMajorAxis();
        std::sort(
                m_tree[nodeIndex].indices.begin(),
                m_tree[nodeIndex].indices.end(),
                [this, dimension](uint32_t f1, uint32_t f2) {
                    return m_mesh->getBoundingBox(f1).getCenter()[dimension] <
                           m_mesh->getBoundingBox(f2).getCenter()[dimension];
                }
        );
        uint32_t index_bucket = 1;
        float cost_min = std::numeric_limits<float>::infinity(); //分桶
        for (uint32_t i = 1; i < COUNT_BUCKET; i++) {
            auto begin_bucket = m_tree[nodeIndex].indices.begin();
            auto mid_bucket = m_tree[nodeIndex].indices.begin() + (static_cast<uint32_t>(m_tree[nodeIndex].indices.size()) * i / COUNT_BUCKET);
            auto end_bucket = m_tree[nodeIndex].indices.end();
            auto faces_left = std::vector<uint32_t>(begin_bucket, mid_bucket);
            auto faces_right = std::vector<uint32_t>(mid_bucket, end_bucket);

            BoundingBox3f bbox_left, bbox_right;
            for (auto left: faces_left) {
                bbox_left = BoundingBox3f::merge(bbox_left, m_mesh->getBoundingBox(left));
            }
            for (auto right: faces_right) {
                bbox_right = BoundingBox3f::merge(bbox_right, m_mesh->getBoundingBox(right));
            }

            auto S_LEFT = bbox_left.getSurfaceArea();
            auto S_RIGHT = bbox_right.getSurfaceArea();
            auto S_ALL = m_tree[nodeIndex].bbox.getSurfaceArea();
            auto cost = 0.125f +
                        static_cast<float>(faces_left.size()) * S_LEFT / S_ALL +
                        static_cast<float>(faces_right.size()) * S_RIGHT / S_ALL;
            if (cost < cost_min) {
                cost_min = cost;
                index_bucket = i;
            }
        }
        auto begin = m_tree[nodeIndex].indices.begin();
        auto mid = m_tree[nodeIndex].indices.begin() + (static_cast<uint32_t>(m_tree[nodeIndex].indices.size()) * index_bucket / COUNT_BUCKET);
        auto end = m_tree[nodeIndex].indices.end();
        auto faces_left = std::vector<uint32_t>(begin, mid);
        auto faces_right = std::vector<uint32_t>(mid, end);

        BoundingBox3f bbox_left, bbox_right;
        for (auto left: faces_left) {
            bbox_left.expandBy(m_mesh->getBoundingBox(left));
        }
        for (auto right: faces_right) {
            bbox_right.expandBy(m_mesh->getBoundingBox(right));
        }

        AccelNode leftNode(bbox_left);
        AccelNode rightNode(bbox_right);
        leftNode.indices = faces_left;
        rightNode.indices = faces_right;

        q.push(m_tree.size());
        m_tree.emplace_back(leftNode);
        q.push(m_tree.size());
        m_tree.emplace_back(rightNode);

        m_count_node += 2;
        m_count_leaf += 1;
    }


    bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
        bool foundIntersection = false;  // Was an intersection found so far?
        auto f = (uint32_t) -1;      // Triangle index of the closest intersection

        Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)


        if (mode == AccelMode::NONE) {
            //Brute force search through all triangles
            for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
                float u, v, t;
                if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
                    //An intersection was found! Can terminate
                    //immediately if this is a shadow ray query
                    if (shadowRay)
                        return true;
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    f = idx;
                    foundIntersection = true;
                }
            }
        } else if (mode == AccelMode::SAH) {
            foundIntersection = traverseBVH(0, ray, its, f, shadowRay);
        } else if (mode == AccelMode::OctTree) {
            foundIntersection = traverseOctTree(0, ray, its, f, shadowRay);
        }


        if (foundIntersection) {
            /* At this point, we now know that there is an intersection,
               and we know the triangle index of the closest such intersection.

               The following computes a number of additional properties which
               characterize the intersection (normals, texture coordinates, etc..)
            */

            /* Find the barycentric coordinates */
            Vector3f bary;
            bary << 1 - its.uv.sum(), its.uv;

            /* References to all relevant mesh buffers */
            const Mesh *mesh = its.mesh;
            const MatrixXf &V = mesh->getVertexPositions();
            const MatrixXf &N = mesh->getVertexNormals();
            const MatrixXf &UV = mesh->getVertexTexCoords();
            const MatrixXu &F = mesh->getIndices();

            /* Vertex indices of the triangle */
            uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

            Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

            /* Compute the intersection positon accurately
               using barycentric coordinates */
            its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

            /* Compute proper texture coordinates if provided by the mesh */
            if (UV.size() > 0)
                its.uv = bary.x() * UV.col(idx0) +
                         bary.y() * UV.col(idx1) +
                         bary.z() * UV.col(idx2);

            /* Compute the geometry frame */
            its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

            if (N.size() > 0) {
                /* Compute the shading frame. Note that for simplicity,
                   the current implementation doesn't attempt to provide
                   tangents that are continuous across the surface. That
                   means that this code will need to be modified to be able
                   use anisotropic BRDFs, which need tangent continuity */

                its.shFrame = Frame(
                        (bary.x() * N.col(idx0) +
                         bary.y() * N.col(idx1) +
                         bary.z() * N.col(idx2)).normalized());
            } else {
                its.shFrame = its.geoFrame;
            }
        }

        return foundIntersection;
    }

    bool Accel::traverseOctTree(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const {
        auto &node = m_tree[n];

        //当前节点包围盒与射线碰撞
        if (!node.bbox.rayIntersect(ray)) {
            return false;
        }

        bool isHit = false;

        //节点为叶子节点
        if (node.child == 0) {
            float u, v, t;
            //遍历节点内的图元
            for (auto i: node.indices) {
                //求与射线相交的最近的图元
                if (m_mesh->rayIntersect(i, ray, u, v, t) && t < ray.maxt) {
                    if (shadowRay) {
                        return true;
                    }
                    ray.maxt = t;
                    its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    f = i;
                    isHit = true;
                }
            }
        } else {
            std::pair<uint32_t, float> children[8] = {};
            for (uint32_t i = 0; i < 8; i++) {
                auto ptr = node.child + i;
                //求出节点到光线原点的距离
                children[i] = {ptr, m_tree[ptr].bbox.distanceTo(ray.o)};
            }
            //按子节点距离排序
            std::sort(
                    children,
                    children + 8,
                    [](const auto &lhs, const auto &rhs) {
                        return lhs.second < rhs.second;
                    }
            );
            //递归子节点的求交
            for (auto &child: children) {
                isHit |= traverseOctTree(child.first, ray, its, f, shadowRay);
                if (shadowRay && isHit) {
                    return true;
                }
            }
        }
        return isHit;
    }

    bool Accel::traverseBVH(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const {
        auto &node = m_tree[n];

        //当前节点包围盒与射线碰撞
        if (!node.bbox.rayIntersect(ray)) {
            return false;
        }

        bool isHit = false;

        //节点为叶子节点
        if (node.child == 0) {
            float u, v, t;
            //遍历节点内的图元
            for (auto i: node.indices) {
                //求与射线相交的最近的图元
                if (m_mesh->rayIntersect(i, ray, u, v, t) && t < ray.maxt) {
                    if (shadowRay) {
                        return true;
                    }
                    ray.maxt = t;
                    its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    f = i;
                    isHit = true;
                }
            }
        } else {
            isHit |= traverseBVH(node.child, ray, its, f, shadowRay);
            isHit |= traverseBVH(node.child + 1, ray, its, f, shadowRay);
            if (shadowRay && isHit) {
                return true;
            }
        }
        return isHit;
    }

NORI_NAMESPACE_END

