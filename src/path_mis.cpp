#include "nori/integrator.h"
#include "nori/scene.h"
#include "nori/bsdf.h"
#include "nori/emitter.h"
#include "nori/sampler.h"

NORI_NAMESPACE_BEGIN

class PathMultiImportanceSampleIntegrator : public Integrator {
public:
    PathMultiImportanceSampleIntegrator(const PropertyList& props) {
        m_maxDepth = props.getInteger("max_depth", 128);
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Color3f L(0.0f);
        Color3f beta(1.0f);
        Intersection its;
        Ray3f ray_(ray);
        int depth = 0;
        bool specularBounce = false;

        while (depth < m_maxDepth)
        {
            // miss
            if (!scene->rayIntersect(ray_, its)) {
                break;
            }

            if (depth == 0 || specularBounce) {
                // hit area light
                L += its.Le(-ray_.d) * beta;
            }
            const BSDF* bsdf = its.mesh->getBSDF();
            if (!bsdf) {
                break;
            }

            BSDFQueryRecord record(its.wi);

            // MIS
            if (bsdf->isDiffuse()) {
                L += beta * uniformSampleOneLight(its, *scene, *sampler);
            }

            // sample by bsdf
            Color3f weight = bsdf->sample(record, sampler->next2D());
            if (weight.isZero()) {
                break;
            }
            beta *= weight;
            specularBounce = record.measure == EDiscrete;

            // russian roulette
            if (depth > 3) {
                float continue_probability = std::min(0.99f, beta.maxCoeff());
                if (sampler->next1D() > continue_probability) {
                    break;
                }
                beta /= continue_probability;
            }
            ray_ = its.SpawnRay(its.toWorld(record.wo));

            ++depth;
        }

        return L;
    }

    std::string toString() const override {
        return "PathMultiImportanceSampleIntegrator []";
    }

private:
    int m_maxDepth = 128;
};

NORI_REGISTER_CLASS(PathMultiImportanceSampleIntegrator, "path_mis");

NORI_NAMESPACE_END