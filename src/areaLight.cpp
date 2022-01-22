#include "nori/emitter.h"
#include "nori/mesh.h"
#include "nori/object.h"
#include "nori/bsdf.h"

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList& props) {
        m_Lemit = props.getColor("radiance");
    }

    Color3f Li(const Intersection& its, const Vector3f& w) const override {
        return its.shFrame.n.dot(w) > 0 ? m_Lemit : Color3f(0.0f);
    }

    Color3f SampleLi(const Intersection& ref, const Point2f& sample, BSDFQueryRecord& bRec, float& pdf, VisibilityTester& vis) const override {
        Intersection its = m_mesh->sample(ref, sample, pdf);
        if (pdf == 0.0f) {
            return Color3f(0.0f);
        }
        Vector3f w = (its.p - ref.p).normalized();
        bRec.wo = ref.toLocal(w);
        bRec.measure = ESolidAngle;
        vis = VisibilityTester(ref, its);
        return this->Li(its, -w);
    }

    float PdfLi(const Intersection& ref, const Vector3f& wo) const override {
        return m_mesh->pdf(ref, wo);
    }

    std::string toString() const override {
        return tfm::format("AreaLight [%s]", m_Lemit.toString());
    }

    void setParent(NoriObject* parent) override {
        if (parent->getClassType() != EMesh) {
            throw NoriException("AreaLight must attch to Mesh, not to %s", parent->toString());
        }
        m_mesh = dynamic_cast<const Mesh*>(parent);
    }

private:
    const Mesh* m_mesh;
    Color3f m_Lemit;
};

NORI_REGISTER_CLASS(AreaLight, "area");

NORI_NAMESPACE_END