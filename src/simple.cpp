#include "nori/integrator.h"
#include "nori/scene.h"

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList& props) {
        m_position = props.getPoint("position");
        m_intensity = props.getColor("energy");
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0.0f);
        }

        Vector3f toLightDirWorld = (m_position - its.p).normalized();
        if (scene->rayIntersect(Ray3f(its.p, toLightDirWorld))) {
            return Color3f(0.0f);
        }

        Vector3f toLightDirLocal = its.shFrame.toLocal(m_position).normalized();
        float cosTheta = std::max(0.0f, its.shFrame.cosTheta(toLightDirLocal));
        float squaredDistance = (m_position - its.p).squaredNorm();

        Color3f luminace = m_intensity * cosTheta / (M_PI * M_PI * 4.0f * squaredDistance);
        return luminace;
    }

    std::string toString() const override {
        return "SimpleIntegrator[]";
    }

private:
    Point3f m_position;
    Color3f m_intensity;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");

NORI_NAMESPACE_END