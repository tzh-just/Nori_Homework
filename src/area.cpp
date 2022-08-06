
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

    class AreaLight : public Emitter
    {
    private:
        Color3f m_radiance;
    public:
        AreaLight(const PropertyList& propList)
        {
            m_radiance = propList.getColor("radiance");
        }

        Color3f eval(const EmitterQueryRecord& record) const override
        {
            return (record.n.dot(record.wi) < 0.0f) ? m_radiance : 0.0f;//法线和wi要面朝不同方向，这样光源正面才朝着着色点
        }

        Color3f getRadiance() const override
        {
            return m_radiance;
        }

        Color3f sample(const Mesh* mesh, EmitterQueryRecord& record, Sampler* sampler) const override
        {
            //对本光源均匀采样
            auto result = mesh->sampleSurfaceUniform(sampler);
            record.p = result.p;
            record.n = result.n;
            record.wi = (record.p - record.ref).normalized();
            record.shadowRay = Ray3f(record.ref, record.wi, Epsilon, (record.p - record.ref).norm() - Epsilon);
            record.pdf = pdf(mesh, record);
            if (record.pdf > 0.0f && !std::isnan(record.pdf) && !std::isinf(record.pdf))
            {
                return eval(record) / record.pdf;
            }
            return Color3f(0.0f);
        }

        float pdf(const Mesh* mesh, const EmitterQueryRecord& record) const override
        {
            float costTheta = record.n.dot(-record.wi);
            if (costTheta <= 0.0f)
            {
                return 0.0f;
            }
            //光源pdf转为立体角上
            return mesh->getPdf().getNormalization() * (record.p - record.ref).squaredNorm() / costTheta;
        }

        std::string toString() const override
        {
            return "Emitter[]";
        }
    };

    NORI_REGISTER_CLASS(AreaLight, "area")
NORI_NAMESPACE_END