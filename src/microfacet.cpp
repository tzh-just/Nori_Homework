/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN
    static float DistributeBeckmann(const Vector3f& wh, float alpha)
    {
        //Beckmann法线分布项
        float tanTheta = Frame::tanTheta(wh);
        float cosTheta = Frame::cosTheta(wh);
        float a = std::exp(-(tanTheta * tanTheta) / (alpha * alpha));
        float b = M_PI * alpha * alpha * std::pow(cosTheta, 4.0f);
        return a / b;
    }

    static float G1(const Vector3f& wv, const Vector3f& wh, float alpha)
    {
        //Beckmann几何遮蔽项
        float c = wv.dot(wh) / Frame::cosTheta(wv);
        if (c <= 0)
        {
            return 0;
        }
        float b = 1.0f / (alpha * Frame::tanTheta(wv));
        return b < 1.6f ? (3.535f * b + 2.181f * b * b) / (1.f + 2.276f * b + 2.577f * b * b) : 1;
    }

    class Microfacet : public BSDF
    {
    public:
        Microfacet(const PropertyList& propList)
        {
            /* RMS surface roughness */
            m_alpha = propList.getFloat("alpha", 0.1f);

            /* Interior IOR (default: BK7 borosilicate optical glass) */
            m_intIOR = propList.getFloat("intIOR", 1.5046f);

            /* Exterior IOR (default: air) */
            m_extIOR = propList.getFloat("extIOR", 1.000277f);

            /* Albedo of the diffuse base material (a.k.a "kd") */
            m_kd = propList.getColor("kd", Color3f(0.5f));

            /* To ensure energy conservation, we must scale the
               specular component by 1-kd.

               While that is not a particularly realistic model of what
               happens in reality, this will greatly simplify the
               implementation. Please see the course staff if you're
               interested in implementing a more realistic version
               of this BRDF. */
            m_ks = 1 - m_kd.maxCoeff();
        }

        /// Evaluate the BRDF for the given pair of directions
        Color3f eval(const BSDFQueryRecord& bRec) const
        {
            Vector3f wh = (bRec.wi + bRec.wo).normalized();
            auto cosThetaI = Frame::cosTheta(bRec.wi);
            auto cosThetaO = Frame::cosTheta(bRec.wo);
            auto D = DistributeBeckmann(wh, m_alpha);
            auto G = G1(bRec.wi, wh, m_alpha)*G1(bRec.wo, wh, m_alpha);
            auto F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
            return m_kd / M_PI + m_ks * ((D * F * G) / (4.0f * cosThetaI * cosThetaO));
        }

        /// Evaluate the sampling density of \ref sample() wrt. solid angles
        float pdf(const BSDFQueryRecord& bRec) const
        {
            if (Frame::cosTheta(bRec.wo) <= 0)
            {
                return 0;
            }

            Vector3f wh = (bRec.wi + bRec.wo).normalized();
            auto D = DistributeBeckmann(wh, m_alpha);
            float jacobian = 1 / (4.0f * abs(wh.dot(bRec.wo)));
            return m_ks * D * Frame::cosTheta(wh) * jacobian + (1 - m_ks) * Frame::cosTheta(bRec.wo) * INV_PI;
        }

        /// Sample the BRDF
        Color3f sample(BSDFQueryRecord& bRec, const Point2f& _sample) const
        {
            if (Frame::cosTheta(bRec.wi) <= 0) {
                return Color3f(0.0f);
            }
            if (_sample.x() > m_ks) {  //漫反射
                Point2f sample((_sample.x() - m_ks) / (1.f - m_ks), _sample.y());
                bRec.wo = Warp::squareToCosineHemisphere(sample);
            } else {  //高光
                Point2f sample(_sample.x() / m_ks, _sample.y());
                Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
                bRec.wo = ((2.0f * wh.dot(bRec.wi) * wh) - bRec.wi).normalized();
            }
            if (Frame::cosTheta(bRec.wo) < 0.f) {
                return Color3f(0.0f);
            }

            return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        }

        bool isDiffuse() const
        {
            /* While microfacet BRDFs are not perfectly diffuse, they can be
               handled by sampling techniques for diffuse/non-specular materials,
               hence we return true here */
            return true;
        }

        std::string toString() const
        {
            return tfm::format(
                    "Microfacet[\n"
                    "  alpha = %f,\n"
                    "  intIOR = %f,\n"
                    "  extIOR = %f,\n"
                    "  kd = %s,\n"
                    "  ks = %f\n"
                    "]",
                    m_alpha,
                    m_intIOR,
                    m_extIOR,
                    m_kd.toString(),
                    m_ks
            );
        }

    private:
        float m_alpha;
        float m_intIOR, m_extIOR;
        float m_ks;
        Color3f m_kd;
    };

    NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
