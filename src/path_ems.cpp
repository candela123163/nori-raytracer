#include "nori/integrator.h"
#include "nori/scene.h"
#include "nori/bsdf.h"
#include "nori/emitter.h"
#include "nori/sampler.h"

NORI_NAMESPACE_BEGIN

class PathEmitterSampleIntegrator : public Integrator {
public:
    PathEmitterSampleIntegrator(const PropertyList& props) {
        m_maxDepth = props.getInteger("max_depth", 128);
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Color3f L(0.0f);
        Color3f beta(1.0f);
        Intersection its;
        Ray3f ray_(ray);
        int depth = 0;
        bool specularBounce = false;

        while(depth < m_maxDepth)
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

            BSDFQueryRecord record(its.toLocal(-ray_.d.normalized()));

            // do direct illumination by sampling light
            if (bsdf->isDiffuse()) {
                // uniformally sample one light
                size_t lightCount = scene->getEmitters().size();
                size_t lightIndex = std::min(lightCount - 1, static_cast<size_t>(sampler->next1D() * lightCount));
                const Emitter* emitter = scene->getEmitters()[lightIndex];
                float lightPdf = 1.0f / lightCount;

                float pdf = 0.0f;
                VisibilityTester vis;
                Color3f Li = emitter->SampleLi(its, sampler->next2D(), record, pdf, vis);
                if (!Li.isZero() && pdf > 0.0f)
                {
                    Color3f weight = bsdf->eval(record) * std::max(0.0f, Frame::cosTheta(record.wo)) / (pdf * lightPdf);
                    if (!weight.isZero() && vis.Unoccluded(*scene)) {
                        L += weight * Li * beta;
                    }
                }
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
                float continue_probility = std::min(0.99f, beta.maxCoeff());
                if (sampler->next1D() > continue_probility) {
                    break;
                }
                beta /= continue_probility;
            }
            ray_ = its.SpawnRay(its.toWorld(record.wo));

            ++depth;
        }

        return L;
    }

    std::string toString() const override {
        return "PathMatSampleIntegrator []";
    }

private:
    int m_maxDepth = 128;
};

NORI_REGISTER_CLASS(PathEmitterSampleIntegrator, "path_ems");

NORI_NAMESPACE_END