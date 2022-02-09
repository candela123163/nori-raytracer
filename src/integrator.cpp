#include "nori/integrator.h"
#include "nori/scene.h"
#include "nori/sampler.h"
#include "nori/emitter.h"
#include "nori/bsdf.h"

NORI_NAMESPACE_BEGIN

float Integrator::powerHeuristic(int nf, float fPdf, int ng, float gPdf) {
	float f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}

float Integrator::balanceHeuristic(int nf, float fPdf, int ng, float gPdf) {
	float f = nf * fPdf, g = ng * gPdf;
	return f / (f + g);
}

Color3f Integrator::uniformSampleOneLight(const Intersection& its, const Scene& scene, Sampler& sampler) {
	size_t lightCount = scene.getEmitters().size();
	size_t lightIndex = std::min(lightCount - 1, static_cast<size_t>(sampler.next1D() * lightCount));
	const Emitter* emitter = scene.getEmitters()[lightIndex];
	float lightPdf = 1.0f / lightCount;
	return estimateDirect(its, scene, sampler, *emitter) / lightPdf;
}

Color3f Integrator::estimateDirect(const Intersection& its, const Scene& scene, Sampler& sampler, const Emitter& emitter) {
	Color3f Ld(0.0f);
	float lightPdf = 0.0f, scatteringPdf = 0.0f;
	VisibilityTester vis;
	BSDFQueryRecord record(its.wi);
	const BSDF* bsdf = its.mesh->getBSDF();

	// sample light
	Color3f Li = emitter.SampleLi(its, sampler.next2D(), record, lightPdf, vis);
	if (lightPdf > 0.0f && !Li.isZero()) {
		Color3f f = bsdf->eval(record) * std::abs(Frame::cosTheta(record.wo));
		scatteringPdf = bsdf->pdf(record);
		if (scatteringPdf > 0.0f && !f.isZero() && vis.Unoccluded(scene)) {
			float weight = powerHeuristic(1, lightPdf, 1, scatteringPdf);
			Ld += Li * f * weight / lightPdf;
		}
	}

	// sample bsdf
	Color3f f = bsdf->sample(record, sampler.next2D());
	scatteringPdf = bsdf->pdf(record);
	if (scatteringPdf > 0.0f && !f.isZero()) {
		Vector3f worldWo = its.toWorld(record.wo);
		lightPdf = emitter.PdfLi(its, worldWo);
		if (lightPdf > 0.0f) {
			Intersection it;
			if (scene.rayIntersect(its.SpawnRay(worldWo), it)) {
				if (it.mesh->getEmitter() == &emitter)
				{
					Color3f Li = it.Le(-worldWo);
					if (!Li.isZero()) {
						float weight = powerHeuristic(1, scatteringPdf, 1, lightPdf);
						Ld += f * Li * weight;
					}
				}
			}
		}
	}

	return Ld;
}

NORI_NAMESPACE_END