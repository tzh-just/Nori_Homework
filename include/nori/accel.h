/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <queue>

#include <nori/mesh.h>

#include <utility>

NORI_NAMESPACE_BEGIN

    enum AccelMode {
        NONE = 0,
        OctTree = 1,
        BVH = 2,
        SAH = 3
    };

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

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
    class Accel {
    public:
        /**
         * \brief Register a triangle mesh for inclusion in the acceleration
         * data structure
         *
         * This function can only be used before \ref build() is called
         */
        void addMesh(Mesh *mesh);

        /// Build the acceleration data structure (currently a no-op)
        void build();

        void buildOctTree(std::queue<uint32_t> &q, uint32_t nodeIndex);

        void buildSAH(std::queue<uint32_t> &q, uint32_t nodeIndex);


        bool traverseOctTree(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const;

        bool traverseBVH(uint32_t n, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const;

        /// Return an axis-aligned box that bounds the scene
        const BoundingBox3f &getBoundingBox() const { return m_bbox; }

        /**
         * \brief Intersect a ray against all triangles stored in the scene and
         * return detailed intersection information
         *
         * \param ray
         *    A 3-dimensional ray data structure with minimum/maximum extent
         *    information
         *
         * \param its
         *    A detailed intersection record, which will be filled by the
         *    intersection query
         *
         * \param shadowRay
         *    \c true if this is a shadow ray query, i.e. a query that only aims to
         *    find out whether the ray is blocked or not without returning detailed
         *    intersection information.
         *
         * \return \c true if an intersection was found
         */
        bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

    private:
        Mesh *m_mesh = nullptr; ///< Mesh (only a single one for now)
        BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
        std::vector<AccelNode> m_tree;

        uint32_t m_depth_tree = 1;
        uint32_t m_count_leaf = 1;
        uint32_t m_count_node = 1;

        static constexpr uint32_t COUNT_MIN = 16;
        static constexpr uint32_t DEPTH_OCT_MAX = 12;
        static constexpr uint32_t DEPTH_BVH_MAX = 32;
        static constexpr uint32_t COUNT_BUCKET = 10;

        AccelMode mode = AccelMode::OctTree;
    };

NORI_NAMESPACE_END
