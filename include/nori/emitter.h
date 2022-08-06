/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

    struct EmitterQueryRecord
    {
        Point3f ref;
        Point3f p;
        Normal3f n;
        Vector3f wi;
        float pdf;
        Ray3f shadowRay;

        EmitterQueryRecord(const Point3f& ref) : ref(ref) {}

        EmitterQueryRecord(const Point3f& ref, const Point3f& p, const Normal3f& n) : ref(ref), p(p), n(n)
        {
            wi = (p - ref).normalized();
        }
    };

/**
 * \brief Superclass of all emitters
 */
    class Emitter : public NoriObject
    {
    public:
        virtual ~Emitter() {}

        virtual Color3f eval(const EmitterQueryRecord& record) const = 0;

        virtual float pdf(const Mesh* mesh, const EmitterQueryRecord& lRec) const = 0;

        virtual Color3f sample(const Mesh* mesh, EmitterQueryRecord& lRec, Sampler*) const = 0;

        virtual Color3f getRadiance() const = 0;

        /**
         * \brief Return the type of object (i.e. Mesh/Emitter/etc.)
         * provided by this instance
         * */
        EClassType getClassType() const { return EEmitter; }
    };

NORI_NAMESPACE_END
