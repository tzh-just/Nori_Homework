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

  virtual void divide(uint32_t n, std::vector<AccelNode> *children) = 0;

  virtual bool traverse(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const = 0;

  virtual std::pair<uint32_t,uint32_t> getLimits() const = 0;

  const BoundingBox3f &getBoundingBox() { return m_bbox; }

 protected:
  Mesh *m_mesh = nullptr;
  BoundingBox3f m_bbox;
  std::vector<AccelNode> m_tree;

  uint32_t depth_curr_ = 1;
  uint32_t count_leaf_ = 1;
  uint32_t count_node_ = 1;
};

class BVH : public Accel {
 public:
  void divide(uint32_t n, std::vector<AccelNode> *children) override;
  bool traverse(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const override;
  std::pair<uint32_t,uint32_t> getLimits() const override;
 private:
  uint32_t COUNT_MIN = 16;
  uint32_t DEPTH_MAX = 32;
  uint32_t COUNT_BUCKET = 10;
};

class OctTree : public Accel {
 public:
  void divide(uint32_t n, std::vector<AccelNode> *children) override;
  bool traverse(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const override;
  std::pair<uint32_t,uint32_t> getLimits() const override;
 private:
  uint32_t COUNT_MIN = 16;
  uint32_t DEPTH_MAX = 12;
};

NORI_NAMESPACE_END
