#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN
    class WhittedIntegrator : public Integrator
    {
    public:
        WhittedIntegrator(const PropertyList& props) {}

        Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
        {
            Intersection its;
            Color3f color(0.0f);
            //判断击中场景中物体
            if(!scene->rayIntersect(ray,its)){
                return color;
            }
            Color3f Le(0.0f);
            //直接击中光源
            if(its.mesh->isEmitter()){
                EmitterQueryRecord record(ray.o,its.p,its.shFrame.n);
                Le = its.mesh->getEmitter()->eval(record);
            }
            //未击中则计算直接光照
            Mesh* light = scene->getRandomEmitter(sampler->next1D());
            EmitterQueryRecord directRecord(its.p);
            Color3f Li = light->getEmitter()->sample(light,directRecord,sampler);
            if(scene->rayIntersect(directRecord.shadowRay)){
                Li = 0;
            }
            float cosTheta = Frame::cosTheta(its.shFrame.toLocal(directRecord.wi));
            BSDFQueryRecord bsdf(its.toLocal(-ray.d),its.toLocal(directRecord.wi),ESolidAngle);
            Color3f f = its.mesh->getBSDF()->eval(bsdf);
            if(cosTheta<0){
                cosTheta = 0;
            }
            return Le+Li*f*cosTheta/(1.0f/(float)scene->getEmitters().size());
        }

        std::string toString() const
        {
            return "WhittedIntegrator[]";
        }
    };

    NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END