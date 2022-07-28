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

std::pair<uint32_t, uint32_t> OctTree::getLimits() const {
  return {16, 12};//count,depth
}

std::pair<uint32_t, uint32_t> BVH::getLimits() const {
  return {16, 32};//count,depth
}

void Accel::build() {
  if (!m_mesh) return;

  auto start = std::chrono::high_resolution_clock::now();

  //初始化树
  m_tree.clear();
  auto root = AccelNode(m_mesh->getBoundingBox(), m_mesh->getTriangleCount());
  for (int i = 0; i < root.indices.size(); i++) {
    root.indices[i] = i;
  }
  m_tree.emplace_back(root);

  //初始化辅助队列
  std::queue<uint32_t> q;
  q.push(0);

  //初始化子节点列表
  auto children = std::vector<AccelNode>();

  //获取限制条件
  auto [LIMIT_COUNT, LIMIT_DEPTH] = getLimits();
  //构建树
  while (!q.empty()) {
    uint32_t n = q.size();//层次遍历
    for (uint32_t k = 0; k < n; k++) {
      //限制之内对节点继续分叉
      if (m_tree[q.front()].indices.size() > LIMIT_COUNT && depth_curr_ < LIMIT_DEPTH) {
        //设置子节点起始索引
        m_tree[q.front()].child = m_tree.size();
        //分割子节点
        divide(q.front(), &children);
        //子节点加入树，索引入队
        --count_leaf_;
        for (auto &child : children) {
          q.push(m_tree.size());
          m_tree.emplace_back(child);
          ++count_node_;
          ++count_leaf_;
        }
      }
      //清理无用数据
      q.pop();
      children.clear();
      children.shrink_to_fit();
    }
    depth_curr_++;
  }

  auto over = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(over - start).count();

  std::cout << "[build time]: " << time << "ms" << std::endl;
  std::cout << "[max depth]: " << depth_curr_ << std::endl;
  std::cout << "[node count]: " << count_node_ << std::endl;
  std::cout << "[leaf count]: " << count_leaf_ << std::endl;
}

void OctTree::divide(uint32_t n, std::vector<AccelNode> *children) {
  auto & node =  m_tree[n];
  //获取包围盒中心点
  Vector3f center = node.bbox.getCenter();
  //获取包围盒八个拐角点
  for (size_t i = 0; i < 8; i++) {
    //构建子包围盒
    Vector3f corner = node.bbox.getCorner(i);
    BoundingBox3f bbox_sub;
    for (uint32_t j = 0; j < 3; j++) {
      bbox_sub.min[j] = std::min(center[j], corner[j]);
      bbox_sub.max[j] = std::max(center[j], corner[j]);
    }

    //构建子节点
    AccelNode node_sub(bbox_sub);
    for (auto face : node.indices) {
      //检测节点持有的图元与子包围盒是否重叠
      if (bbox_sub.overlaps(m_mesh->getBoundingBox(face))) {
        node_sub.indices.emplace_back(face);
      }
    }
    children->emplace_back(node_sub);
  }
}

void BVH::divide(uint32_t n, std::vector<AccelNode> *children) {
  auto & node =  m_tree[n];
  //在最长维度排序
  uint32_t dimension = node.bbox.getMajorAxis();
  std::sort(
      node.indices.begin(),
      node.indices.end(),
      [this, dimension](uint32_t f1, uint32_t f2) {
        return m_mesh->getBoundingBox(f1).getCenter()[dimension] <
            m_mesh->getBoundingBox(f2).getCenter()[dimension];
      }
  );

  //分桶
  float cost_min = std::numeric_limits<float>::infinity();
  AccelNode leftNode, rightNode;
  std::vector<uint32_t> faces_left, faces_right;
  for (uint32_t i = 1; i < COUNT_BUCKET; i++) {
    //根据通的数量依次划分左右
    auto begin = node.indices.begin();
    auto mid = node.indices.begin() + (static_cast<uint32_t>(node.indices.size()) * i / COUNT_BUCKET);
    auto end = node.indices.end();
    faces_left = std::vector<uint32_t>(begin, mid);
    faces_right = std::vector<uint32_t>(mid, end);

    //合并包围盒
    BoundingBox3f bbox_left, bbox_right;
    for (auto left : faces_left) {
      bbox_left.expandBy(m_mesh->getBoundingBox(left));
    }
    for (auto right : faces_right) {
      bbox_right.expandBy(m_mesh->getBoundingBox(right));
    }

    //计算成本
    float S_LEFT = bbox_left.getSurfaceArea();
    float S_RIGHT = bbox_right.getSurfaceArea();
    float S_ALL = node.bbox.getSurfaceArea();
    float cost = 0.125f +
        static_cast<float>(faces_left.size()) * S_LEFT / S_ALL +
        static_cast<float>(faces_right.size()) * S_RIGHT / S_ALL;

    //根据成本选择分桶方案
    if (cost < cost_min) {
      cost_min = cost;
      leftNode.bbox = std::move(bbox_left);
      leftNode.indices = std::move(faces_left);
      rightNode.bbox = std::move(bbox_right);
      rightNode.indices = std::move(faces_right);
    }
  }
  children->emplace_back(leftNode);
  children->emplace_back(rightNode);
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
  bool foundIntersection = false;  // Was an intersection found so far?
  auto f = (uint32_t) -1;      // Triangle index of the closest intersection

  Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)


  //Brute force search through all triangles
/*  for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
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
  }*/

  foundIntersection = traverse(0, ray, its, f, shadowRay);

  if (shadowRay)
    return foundIntersection;

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

bool OctTree::traverse(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const {
  auto &node = m_tree[n];

  //当前节点包围盒与射线碰撞
  if (!node.bbox.rayIntersect(ray))
    return false;

  bool isHit = false;

  //节点为叶子节点
  if (node.child == 0) {
    float u, v, t;
    //遍历节点内的图元
    for (auto i : node.indices) {
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
    for (auto &child : children) {
      isHit |= traverse(child.first, ray, its, f, shadowRay);
      if (shadowRay && isHit) {
        return true;
      }
    }
  }
  return isHit;
}

bool BVH::traverse(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const {
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
    for (auto i : node.indices) {
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
    isHit |= traverse(node.child, ray, its, f, shadowRay);
    isHit |= traverse(node.child + 1, ray, its, f, shadowRay);
    if (shadowRay && isHit) {
      return true;
    }
  }
  return isHit;
}

NORI_NAMESPACE_END

