/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

//
    static Vector3f refract(const Vector3f& viewDir, const Vector3f& normal, float eta)
    {
        //判断光线在物体内部还是外部
        bool isOuter = normal.dot(viewDir) > 0;
        auto N = isOuter ? normal : -normal;
        eta = isOuter ? eta : (1.0f / eta);

        float cosThetaI = viewDir.dot(N);
        float cosThetaT2 = 1 - eta * eta * (1 - cosThetaI * cosThetaI);

        if (cosThetaT2 <= 0)
        {
            return Vector3f(0.0f);
        }
        float cosThetaT = sqrt(cosThetaT2);
        return -viewDir * eta + N * (cosThetaI * eta - cosThetaT);
    }

/// Ideal dielectric BSDF
    class Dielectric : public BSDF
    {
    public:
        Dielectric(const PropertyList& propList)
        {
            /* Interior IOR (default: BK7 borosilicate optical glass) */
            m_intIOR = propList.getFloat("intIOR", 1.5046f);

            /* Exterior IOR (default: air) */
            m_extIOR = propList.getFloat("extIOR", 1.000277f);
        }

        Color3f eval(const BSDFQueryRecord&) const
        {
            /* Discrete BRDFs always evaluate to zero in Nori */
            return Color3f(0.0f);
        }

        float pdf(const BSDFQueryRecord&) const
        {
            /* Discrete BRDFs always evaluate to zero in Nori */
            return 0.0f;
        }

        Color3f sample(BSDFQueryRecord& bRec, const Point2f& sample) const
        {//throw NoriException("Unimplemented!");
            Vector3f normalLocal{0., 0., 1.};
            float cosThetaI = Frame::cosTheta(bRec.wi);
            bool entering = (cosThetaI > 0.f);
            float etaI, etaT;
            if (entering)
            {
                etaI = m_extIOR;
                etaT = m_intIOR;
            }
            else
            {
                etaI = m_intIOR;
                etaT = m_extIOR;
                normalLocal = -normalLocal;
                cosThetaI = -cosThetaI;
            }
            bRec.eta = etaI / etaT;
            float sin2ThetaI = Frame::sinTheta2(bRec.wi);
            float sin2ThetaT = bRec.eta * bRec.eta * sin2ThetaI;
            if (sin2ThetaT >= 1.)
            { //全反射
                bRec.wo = Vector3f{-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z()};
            }
            else
            {
                //反射的比例，精确
                float fresnelTerm = nori::fresnel(cosThetaI, etaI, etaT);
                //反射的比例，近似
                //float fre =  SchlicksFresnel(cosThetaI, etaI, etaT);
                if (sample.x() < fresnelTerm)
                {
                    bRec.wo = Vector3f{-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z()};
                }
                else
                {
                    float cosThetaT = std::sqrt(1. - sin2ThetaT);
                    //球坐标方法

                    float thetaT;
                    if (entering)
                    {
                        thetaT = M_PI - acos(cosThetaT);
                    }
                    else
                    {
                        thetaT = acos(cosThetaT);
                    }
                    float phiT = std::atan2(bRec.wi.y(), bRec.wi.x()) + M_PI;
                    bRec.wo = sphericalDirection(thetaT, phiT);

                    //向量方法
                    //bRec.wo = bRec.eta * -bRec.wi + (bRec.eta * cosThetaI - cosThetaT) * normalLocal;
                }
            }
            return Color3f(1.);
        }

        std::string toString() const
        {
            return tfm::format(
                    "Dielectric[\n"
                    "  intIOR = %f,\n"
                    "  extIOR = %f\n"
                    "]",
                    m_intIOR, m_extIOR);
        }

    private:
        float m_intIOR, m_extIOR;
    };

    NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
