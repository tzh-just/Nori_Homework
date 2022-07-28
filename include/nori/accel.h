/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <queue>

#include <nori/mesh.h>

#include <utility>

NORI_NAMESPACE_BEGIN
struct AccelNode {
  uint32_t child = 0;
  BoundingBox3f bbox;
  std::vector<uint32_t> indices;

  AccelNode() : bbox() {}
  explicit AccelNode(BoundingBox3f box)
      : bbox(std::move(box)) {}
  AccelNode(BoundingBox3f box, uint32_t size)
      : bbox(std::move(box)), indices(size) {}
};

class Accel {
 public:
  void addMesh(Mesh *mesh);

  void build();

  bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

  virtual bool divide(uint32_t nodeIndex, std::vector<AccelNode> *children) = 0;

  virtual bool traverse(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const = 0;

  const BoundingBox3f &getBoundingBox() { return m_bbox; }

 protected:
  Mesh *m_mesh = nullptr;
  BoundingBox3f m_bbox;
  std::vector<AccelNode> m_tree;

  uint32_t m_depth_tree = 1;
  uint32_t m_count_leaf = 1;
  uint32_t m_count_node = 1;
};

class BVH : public Accel {
 public:
  bool divide(uint32_t nodeIndex, std::vector<AccelNode> *children) override;
  bool traverse(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const override;
 private:
  static constexpr uint32_t COUNT_BVH_MIN = 16;
  static constexpr uint32_t DEPTH_BVH_MAX = 32;
  static constexpr uint32_t COUNT_BUCKET = 10;
};

class OctTree : public Accel {
 public:
  bool divide(uint32_t nodeIndex, std::vector<AccelNode> *children) override;
  bool traverse(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const override;
 private:
  static constexpr uint32_t COUNT_OCT_MIN = 16;
  static constexpr uint32_t DEPTH_OCT_MAX = 12;
};

NORI_NAMESPACE_END
