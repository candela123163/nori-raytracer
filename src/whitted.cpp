#include "nori/integrator.h"
#include "nori/scene.h"
#include "nori/bsdf.h"
#include "nori/emitter.h"
#include "nori/sampler.h"

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList& props) {
        m_maxDepth = props.getInteger("max_depth", 16);
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Color3f L(0.0f);
        Color3f f = 1.0f;
        Intersection its;
        Ray3f ray_(ray);
        int depth = 0;

        while (depth++ < m_maxDepth)
        {
            // miss
            if (!scene->rayIntersect(ray_, its)) {
                break;
            }
            // hit area light
            L += its.Le(-ray_.d);
            const BSDF* bsdf = its.mesh->getBSDF();
            if (!bsdf) {
                break;
            }

            BSDFQueryRecord record(its.toLocal(-ray_.d.normalized()));
            // diffuse
            if (bsdf->isDiffuse()) {
                // sample light
                for (auto emitter : scene->getEmitters())
                {
                    float pdf = 0.0f;
                    VisibilityTester vis;
                    Color3f Li = emitter->SampleLi(its, sampler->next2D(), record, pdf, vis);
                    if (Li.isZero() || pdf == 0.0f) {
                        continue;
                    }
                    Color3f weight = bsdf->eval(record);
                    if (!weight.isZero() && vis.Unoccluded(*scene)) {
                        L += weight * Li * Frame::cosTheta(record.wo) / pdf;
                    }
                }
                break;
            }
            // specular reflection and transmission
            else {
                Color3f weight = bsdf->sample(record, sampler->next2D());
                if (!weight.isZero()) {
                    f *= weight;
                    ray_ = its.SpawnRay(its.toWorld(record.wo));
                }
            }
        }

        return L * f;
    }

    std::string toString() const override {
        return "WhittedIntegrator []";
    }

private:
    int m_maxDepth = 16;
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");

NORI_NAMESPACE_END