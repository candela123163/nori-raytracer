#include "nori/integrator.h"
#include "nori/scene.h"
#include "nori/bsdf.h"
#include "nori/emitter.h"
#include "nori/sampler.h"

NORI_NAMESPACE_BEGIN

class PathMatSampleIntegrator : public Integrator {
public:
    PathMatSampleIntegrator(const PropertyList& props) {
        m_maxDepth = props.getInteger("max_depth", 128);
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Color3f L(0.0f);
        Color3f beta(1.0f);
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
            L += its.Le(-ray_.d) * beta;
            const BSDF* bsdf = its.mesh->getBSDF();
            if (!bsdf) {
                break;
            }

            BSDFQueryRecord record(its.toLocal(-ray_.d.normalized()));
  
            // sample by material's bsdf
            Color3f weight = bsdf->sample(record, sampler->next2D());
            if (weight.isZero()) {
                break;
            }
            beta *= weight;

            // russian roulette
            if (depth > 3) {
                float continue_probility = std::min(0.99f, beta.maxCoeff());
                if (sampler->next1D() > continue_probility) {
                    break;
                }
                beta /= continue_probility;
            }
            ray_ = its.SpawnRay(its.toWorld(record.wo));
        }

        return L;
    }

    std::string toString() const override {
        return "PathMatSampleIntegrator []";
    }

private:
    int m_maxDepth = 128;
};

NORI_REGISTER_CLASS(PathMatSampleIntegrator, "path_mats");

NORI_NAMESPACE_END