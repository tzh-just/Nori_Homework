/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <chrono>
#include <queue>
#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

    void Accel::addMesh(Mesh *mesh) {
        if (m_mesh)
            throw NoriException("Accel: only a single mesh is supported!");
        m_mesh = mesh;
        m_bbox = m_mesh->getBoundingBox();
    }

#define OctTree

    void Accel::build() {
#ifdef BVH
        buildBVH();
#endif

#ifdef OctTree
        buildOctTree();
#endif
    }

    void Accel::buildOctTree() {
        if (!m_mesh) return;

        auto start = std::chrono::high_resolution_clock::now();

        m_octTree.clear();

        auto root = OctNode(m_mesh->getBoundingBox(), m_mesh->getTriangleCount());

        for (int i = 0; i < root.indices.size(); i++) {
            root.indices[i] = i;
        }

        //添加根节点
        m_octTree.emplace_back(root);

        //辅助队列，存储八叉树节点的索引
        std::queue<uint32_t> q;
        q.push(0);

        //记录当前深度
        uint32_t depth_curr = 1;

        //构建八叉树
        while (!q.empty()) {
            //层次遍历
            uint32_t n = q.size();
            for (uint32_t k = 0; k < n; k++) {

                //出队
                auto ptr = q.front();
                q.pop();

                //节点持有的图元数量不足以继续分裂
                if (m_octTree[ptr].indices.size() < COUNT_MIN) {
                    continue;
                }

                //八叉树构建深度超过了最大深度
                if (depth_curr >= DEPTH_MAX) {
                    continue;
                }

                m_octTree[ptr].child = m_octTree.size();


                Vector3f center = m_octTree[ptr].bbox.getCenter();

                //计算节点分裂的八个子节点包围盒
                for (uint32_t i = 0; i < 8; i++) {
                    //获取包围盒中心点和拐角点左边

                    Vector3f corner = m_octTree[ptr].bbox.getCorner(i);

                    //计算子包围盒范围
                    Vector3f minPoint, maxPoint;
                    for (uint32_t j = 0; j < 3; j++) {
                        minPoint[j] = std::min(center[j], corner[j]);
                        maxPoint[j] = std::max(center[j], corner[j]);
                    }

                    //生成子包围盒
                    BoundingBox3f bbox_sub(minPoint, maxPoint);

                    OctNode node_sub(bbox_sub);

                    //遍历节点持有的所有图元
                    for (auto face: m_octTree[ptr].indices)
                        //判断子节点是否覆盖图元
                        if (bbox_sub.overlaps(m_mesh->getBoundingBox(face)))
                            node_sub.indices.emplace_back(face);

                    q.push(m_octTree.size());
                    m_octTree.emplace_back(node_sub);
                }

                m_nodeCount += 8;
                m_leafCount += 7;

                //清理父节点的图元索引数组并销毁无用内存空间
/*                m_octTree[ptr].indices.clear();
                m_octTree[ptr].indices.shrink_to_fit();*/
            }
            depth_curr++;
        }

        m_maxDepth = depth_curr;

        auto end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout << "build time: " << time << std::endl;
        std::cout << "max depth: " << m_maxDepth << std::endl;
        std::cout << "node count: " << m_nodeCount << std::endl;
        std::cout << "leaf count: " << m_leafCount << std::endl;
    }

    void Accel::addOctTreeNode(uint32_t index) {

    }

    void Accel::buildBVH() {
        if (!m_mesh) return;

        auto start = std::chrono::high_resolution_clock::now();

        m_BVH.clear();
        auto root = BVHNode(m_mesh->getBoundingBox());
        m_BVH.emplace_back(root);

        for (uint32_t index = 0; index < m_mesh->getTriangleCount(); index++)
            addBVHNode(index);

        auto end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "build time: " << time << std::endl;
        std::cout << "max depth: " << m_maxDepth << std::endl;
        std::cout << "node count: " << m_nodeCount << std::endl;
        std::cout << "leaf count: " << m_leafCount << std::endl;
    }


    void Accel::addBVHNode(uint32_t index) {

        float cost_min = std::numeric_limits<float>::infinity();

        uint32_t index_min_cost = 0;

        uint32_t bucket = 10;

        //图元的包围盒
        auto box = m_mesh->getBoundingBox(index);
        //辅助队列，存储二叉树节点的索引
        std::queue<uint32_t> q;

        //添加根节点
        q.push(0);

        uint32_t depth = 0;

        //构建二叉树
        while (!q.empty()) {
            uint32_t n = q.size();
            //层次遍历
            for (uint32_t k = 0; k < n; k++) {
                //出队列
                auto ptr = q.front();
                q.pop();

                //节点包围盒未覆盖图元
                if (!m_BVH[ptr].bbox.overlaps(box)) {
                    continue;
                }

                m_maxDepth = std::max(m_maxDepth, depth);

                //节点为叶子节点
                if (m_BVH[ptr].child == 0) {
                    m_BVH[ptr].indices.emplace_back(index);

                    //节点持有的图元数量不足以继续分裂
                    if (m_BVH[ptr].indices.size() < COUNT_MIN) {
                        continue;
                    }

                    //二叉树构建深度超过了最大深度
                    if (depth >= DEPTH_MAX) {
                        continue;
                    }

                    //叶子节点可以作为父节点，并分裂出两个

                    m_nodeCount += 2;
                    m_leafCount += 1;

                    //构造子节点

                    //由相对索引得到节点的绝对索引
                    m_BVH[ptr].child = m_BVH.size();

                    //获取包围盒的最长轴
                    uint32_t dimension = m_BVH[ptr].bbox.getMajorAxis();

                    //按最长维度进行排序
                    std::vector<uint32_t> faces(m_mesh->getTriangleCount());

                    for (int i = 0; i < faces.size(); i++) {
                        faces[i] = i;
                    }

                    switch (dimension) {
                        case 0:
                            std::sort(
                                    faces.begin(),
                                    faces.end(),
                                    [this](uint32_t f1, uint32_t f2) {
                                        return m_mesh->getBoundingBox(f1).getCenter()[0] < m_mesh->getBoundingBox(f2).getCenter()[0];
                                    }
                            );
                            break;
                        case 1:
                            std::sort(
                                    faces.begin(),
                                    faces.end(),
                                    [this](uint32_t f1, uint32_t f2) {
                                        return m_mesh->getBoundingBox(f1).getCenter()[1] < m_mesh->getBoundingBox(f2).getCenter()[1];
                                    }
                            );
                            break;
                        case 2:
                            std::sort(
                                    faces.begin(),
                                    faces.end(),
                                    [this](uint32_t f1, uint32_t f2) {
                                        return m_mesh->getBoundingBox(f1).getCenter()[2] < m_mesh->getBoundingBox(f2).getCenter()[2];
                                    }
                            );
                        default:
                            break;

                    }
                    //在维度上分桶

                    for (uint32_t i = 1; i < bucket; i++) {
                        auto start = faces.begin();
                        auto mid = faces.begin() + faces.size() * i / bucket;
                        auto end = faces.end();


                        auto lefts = std::vector<uint32_t>(start, mid);
                        auto rights = std::vector<uint32_t>(mid, end);


                        BoundingBox3f bbox_left;
                        BoundingBox3f bbox_right;


                        for (auto left: lefts) {
                            bbox_left = BoundingBox3f::merge(bbox_left, m_mesh->getBoundingBox(left));
                        }
                        for (auto right: rights) {
                            bbox_right = BoundingBox3f::merge(bbox_right, m_mesh->getBoundingBox(right));
                        }


                        auto S_left = bbox_left.getSurfaceArea();
                        auto S_right = bbox_left.getSurfaceArea();
                        auto S_all = m_BVH[ptr].bbox.getSurfaceArea();

                        auto cost = 0.125f + lefts.size() * S_left / S_all + rights.size() * S_right / S_all;


                        if (cost < cost_min) {
                            cost_min = cost;
                            index_min_cost = i;
                        }
                    }

                    auto start = faces.begin();
                    auto mid = faces.begin() + faces.size() * index_min_cost / bucket;
                    auto end = faces.end();

                    auto lefts = std::vector<uint32_t>(start, mid);
                    auto rights = std::vector<uint32_t>(mid, end);

                    BoundingBox3f bbox_left;
                    BoundingBox3f bbox_right;

                    for (auto f: lefts) {
                        bbox_left = BoundingBox3f::merge(bbox_left, m_mesh->getBoundingBox(f));
                    }
                    for (auto f: rights) {
                        bbox_right = BoundingBox3f::merge(bbox_right, m_mesh->getBoundingBox(f));
                    }

                    std::cout << bbox_left.toString() << "\n";
                    std::cout << bbox_right.toString() << "\n";

                    BVHNode leftNode(bbox_left);
                    BVHNode rightNode(bbox_right);

                    //图元转移到子节点
                    leftNode.indices = lefts;
                    rightNode.indices = rights;

                    //清理父节点的图元所占内存空间
                    m_BVH[ptr].indices.clear();
                    m_BVH[ptr].indices.shrink_to_fit();

                    m_BVH.emplace_back(leftNode);
                    m_BVH.emplace_back(rightNode);


                } else {
                    //子节点索引都入队列
                    q.push(m_BVH[ptr].child);
                    q.push(m_BVH[ptr].child + 1);
                }
            }
            depth++;
        }
    }

    bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
        bool foundIntersection = false;  // Was an intersection found so far?
        auto f = (uint32_t) -1;      // Triangle index of the closest intersection

        Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

        /* Brute force search through all triangles */
/*        for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
            float u, v, t;
            if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
                *//* An intersection was found! Can terminate
                   immediately if this is a shadow ray query *//*
                if (shadowRay)
                    return true;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                f = idx;
                foundIntersection = true;
            }
        }*/
#ifdef BVH
        foundIntersection = traverseBVH(0, ray, its, f, shadowRay);
#endif

#ifdef OctTree
        foundIntersection = traverseOctTree(0, ray, its, f, shadowRay);
#endif
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
        auto &node = m_octTree[n];

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
                children[i] = {ptr, m_octTree[ptr].bbox.distanceTo(ray.o)};
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
        auto &node = m_BVH[n];

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

