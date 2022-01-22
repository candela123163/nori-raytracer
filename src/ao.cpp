#include "nori/integrator.h"
#include "nori/scene.h"
#include "nori/warp.h"

NORI_NAMESPACE_BEGIN

class AOIntegrator : public Integrator {
public:
    AOIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Color3f L(0.0f);
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return L;
        }

        Point3f shadePoint = its.p;
        Frame shFrame = its.shFrame;
        
        Point3f samplePoint = Warp::squareToCosineHemisphere(sampler->next2D());
        Vector3f sampleDirWorld = shFrame.toWorld(samplePoint);
        if (scene->rayIntersect(Ray3f(shadePoint, sampleDirWorld))) {
            return L;
        }
        float pdf = Warp::squareToCosineHemispherePdf(samplePoint);
        if (pdf < Epsilon) {
            return L;
        }
        float cosTheta = shFrame.cosTheta(samplePoint);
        float f = std::max(0.0f, cosTheta) / (M_PI * pdf);
        L = Color3f(f);
        return L;
    }

    std::string toString() const override {
        return "AOIntegrator []";
    }

};

NORI_REGISTER_CLASS(AOIntegrator, "ao");

NORI_NAMESPACE_END