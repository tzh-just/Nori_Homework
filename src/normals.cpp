#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

    class NormalIntegrator : public Integrator {
    public:
        NormalIntegrator(const PropertyList &props) { }
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
            Intersection its;
            if(!scene->rayIntersect(ray,its)){
                return Color3f(0.0f);
            }
            Normal3f n = its.shFrame.n.cwiseAbs();
            return  Color3f(n.x(),n.y(),n.z());
        }

        std::string toString() const {
            return "NormalIntegrator[]";
        }
    };

    NORI_REGISTER_CLASS(NormalIntegrator, "normals");

NORI_NAMESPACE_END