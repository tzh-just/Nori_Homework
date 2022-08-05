/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

    Point2f Warp::squareToUniformSquare(const Point2f &sample) {
        return sample;
    }

    float Warp::squareToUniformSquarePdf(const Point2f &sample) {
        return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
    }


    Point2f Warp::squareToTent(const Point2f &sample) {
        auto x = sample.x() < 0.5f ? sqrt(2 * sample.x()) - 1 : 1 - sqrt(2 - 2 * sample.x());
        auto y = sample.y() < 0.5f ? sqrt(2 * sample.y()) - 1 : 1 - sqrt(2 - 2 * sample.y());
        return {x, y};
    }

    float Warp::squareToTentPdf(const Point2f &p) {
        bool x = p.x() >= -1 && p.x() <= 1;
        bool y = p.y() >= -1 && p.y() <= 1;

        return x && y ? (1 - abs(p.x())) * (1 - abs(p.y())) : 0;
    }

    Point2f Warp::squareToUniformDisk(const Point2f &sample) {
        auto radius = sqrt(sample.x());
        auto angle = sample.y() * M_PI * 2;
        return {radius * cos(angle), radius * sin(angle)};
    }

    float Warp::squareToUniformDiskPdf(const Point2f &p) {
        return sqrt(p.x() * p.x() + p.y() * p.y()) <= 1 ? INV_PI : 0.0f;
    }

    Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
        auto phi = sample.x() * M_PI * 2;
        auto theta = acos(1 - 2 * sample.y());
        return {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
    }

    float Warp::squareToUniformSpherePdf(const Vector3f &v) {
        return INV_FOURPI;
    }

    Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
        auto phi = sample.x() * M_PI * 2;
        auto theta = acos(1 - sample.y());
        return {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
    }

    float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
        return v.z() >= 0 ? INV_TWOPI : 0;
    }

    Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
        auto d = squareToUniformDisk(sample);
        return {d.x(), d.y(), sqrt(1 - d.x() * d.x() - d.y() * d.y())};
    }

    float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
        return v.z() >= 0 ? v.z()* INV_PI : 0;
    }

    Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
        float phi = M_PI * 2 * sample.x();
        float theta = atan(sqrt(-alpha * alpha * log(1 - sample.y())));
        return {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
    }

    float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
        if (m.z() <= 0)
            return 0;
        float alpha2 = alpha * alpha;
        float tanTheta2 = (m.x() * m.x() + m.y() * m.y()) / (m.z() * m.z());
        float cosTheta3 = m.z() * m.z() * m.z();
        float longitudinal = exp(-tanTheta2 / alpha2) / (alpha2 * cosTheta3);
        return longitudinal * INV_PI;
    }

NORI_NAMESPACE_END
