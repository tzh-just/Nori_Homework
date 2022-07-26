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

            if (scene->rayIntersect(Ray3f(its.p + L * Epsilon, L)))//可见性
                return {0.0f};

            //直接光照：Phi/4pi*pi * cos_theta/||x-p||^2 * V(x<->p)
            return 0.25f * INV_PI * INV_PI * m_energy * std::max(0.0f, its.shFrame.n.dot(L.normalized())) / L.dot(L);
        }

        std::string toString() const {
            return tfm::format(
                    "SimpleIntegrator[\n"
                    "  position = %s\n"
                    "  energy = %s\n"
                    "]", m_position.toString(), m_energy.toString());
        }
    };

    NORI_REGISTER_CLASS(SimpleIntegrator, "simple");

NORI_NAMESPACE_END