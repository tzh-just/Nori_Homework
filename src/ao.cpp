#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

    class AOIntegrator : public Integrator {
    public:
        AOIntegrator(const PropertyList &propList) {}
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
            Intersection its;//相交
            if (!scene->rayIntersect(ray, its))
                return {0.0f};

            auto dir = Warp::squareToCosineHemisphere(sampler->next2D());//cosine-weighted半球采样方向
            auto pdf = Warp::squareToCosineHemispherePdf(dir);
            dir = its.shFrame.toWorld(dir);

            Intersection shadowIts;
            if (scene->rayIntersect(Ray3f(its.p + dir * Epsilon, dir),shadowIts))//可见性
                return {0.0f};

            //Li/pi * cos_theta / pdf
            return Color3f(1.0f) * INV_PI * std::max(0.0f, its.shFrame.n.dot(dir.normalized())) / pdf;
        }

        std::string toString() const {
            return "AOIntegrator[]";
        }
    };
    NORI_REGISTER_CLASS(AOIntegrator, "ao");

NORI_NAMESPACE_END