#include "nori/integrator.h"
#include "nori/scene.h"
#include "nori/warp.h"
#include "pcg32.h"

NORI_NAMESPACE_BEGIN

class AOIntegrator : public Integrator {
public:
    AOIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0.0f);
        }

        pcg32 rng;
        float integration = 0.0f;
        Point3f shadePoint = its.p;
        Frame shFrame = its.shFrame;
        for (uint8_t i = 0; i < m_sampleCount; i++)
        {
            Point2f sample(rng.nextFloat(), rng.nextFloat());
            Point3f samplePoint = Warp::squareToCosineHemisphere(sample);
            Vector3f sampleDirWorld = shFrame.toWorld(samplePoint);
            if (scene->rayIntersect(Ray3f(shadePoint, sampleDirWorld, 0.01f, 3.0f))) {
                continue;
            }
            float pdf = Warp::squareToCosineHemispherePdf(samplePoint);
            if (pdf < Epsilon) {
                continue;
            }
            float cosTheta = shFrame.cosTheta(samplePoint);
            integration += std::max(0.0f, cosTheta) / pdf;
        }
        integration /= m_sampleCount * M_PI;
        return Color3f(integration);
    }

    std::string toString() const override {
        return "AOIntegrator []";
    }

private:
    uint8_t m_sampleCount = 64;
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");

NORI_NAMESPACE_END