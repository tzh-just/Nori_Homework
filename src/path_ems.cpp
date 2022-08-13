#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN
    class PathEmsIntegrator : public Integrator
    {
    private:
        int maxDepth = 10;
    public:
        PathEmsIntegrator(const PropertyList& props) {}

        Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray_) const
        {
            int depth = 0;
            Color3f radiance(0.0f);
            float coff = 1.0f;
            Ray3f ray = ray_;
            while (depth < maxDepth)
            {
                depth++;
                Intersection its;
                //判断击中场景中物体
                if (!scene->rayIntersect(ray, its))
                {
                    continue;
                }
                Color3f Le(0.0f);
                //直接击中光源
                if (its.mesh->isEmitter())
                {
                    EmitterQueryRecord record(ray.o, its.p, its.shFrame.n);
                    Le = its.mesh->getEmitter()->eval(record);
                }
                //漫反射不再弹射
                if (its.mesh->getBSDF()->isDiffuse())
                {
                    //均匀选取光源上一点
                    Mesh* light = scene->getRandomEmitter(sampler->next1D());
                    EmitterQueryRecord directRecord(its.p);
                    //对光源采样
                    Color3f Li = light->getEmitter()->sample(light, directRecord, sampler);

                    //阴影测试
                    if (scene->rayIntersect(directRecord.shadowRay))
                    {
                        Li = 0;
                    }

                    //在着色点采样BxDF
                    float cosTheta = Frame::cosTheta(its.shFrame.toLocal(directRecord.wi));
                    BSDFQueryRecord bsdf(its.toLocal(-ray.d), its.toLocal(directRecord.wi), ESolidAngle);
                    Color3f f = its.mesh->getBSDF()->eval(bsdf);
                    if (cosTheta < 0)
                    {
                        cosTheta = 0;
                    }
                    radiance += Le + Li * f * cosTheta / (1.0f / (float) scene->getEmitters().size());
                    break;
                }
                else
                {
                    //在着色点采样BxDF
                    BSDFQueryRecord bRec(its.toLocal(-ray.d));
                    Color3f refColor = its.mesh->getBSDF()->sample(bRec, sampler->next2D());//采样一个出射方向

                    //俄罗斯轮盘赌
                    if (sampler->next1D() > 0.95 && refColor.x() <= 0.f)
                    {
                        break;
                    }

                    ray = Ray3f(its.p, its.toWorld(bRec.wo));
                    coff /= 0.95f;
                }
            }
            return radiance *= coff;
        }

        std::string toString() const
        {
            return "WhittedIntegrator[]";
        }
    };

    NORI_REGISTER_CLASS(PathEmsIntegrator, "path_ems");
NORI_NAMESPACE_END