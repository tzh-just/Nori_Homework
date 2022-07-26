#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

    class SimpleIntegrator : public Integrator {
    private:
        Point3f m_position;
        Color3f m_energy;
    public:

        SimpleIntegrator(const PropertyList &propList) {
            m_position = propList.getPoint("position", Point3f());
            m_energy = propList.getColor("energy", Color3f());
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
            Intersection its;//相交
            if (!scene->rayIntersect(ray, its))
                return {0.0f};

            Vector3f L = m_position - its.p;//光源方向
            Intersection shadowIts;//可见性
            if (scene->rayIntersect(Ray3f(its.p + L * Epsilon, L), shadowIts))
                return {0.0f};

            //Phi/4pi*pi * max(0, cos_theta)/||x-p||^2 * V(x<->p)
            return 0.25f * INV_PI * INV_PI * m_energy * std::max(0.0f, its.shFrame.n.dot(L.normalized())) / L.dot(L);
        }

        std::string toString() const {
            return tfm::format(
                    "Diffuse[\n"
                    "  position = %s\n"
                    "  energy = %s\n"
                    "]", m_position.toString(), m_energy.toString());
        }
    };

    NORI_REGISTER_CLASS(SimpleIntegrator, "simple");

NORI_NAMESPACE_END