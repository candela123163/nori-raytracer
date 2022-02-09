/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class RoughDieletric : public BSDF {
public:
    RoughDieletric(const PropertyList& propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);
        m_alpha2 = m_alpha * m_alpha;

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the
           specular component by 1-kd.

           While that is not a particularly realistic model of what
           happens in reality, this will greatly simplify the
           implementation. Please see the course staff if you're
           interested in implementing a more realistic version
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    float G1(const Vector3f& wv, const Vector3f& wh) const {
        if (Frame::cosTheta(wv) * wv.dot(wh) <= 0) {
            return 0.0f;
        }
        float b = 1.0f / (m_alpha * std::abs(Frame::tanTheta(wv)));
        if (b >= 1.6) {
            return 1.0f;
        }
        return (3.535 * b + 2.181 * b * b) / (1 + 2.276 * b + 2.577 * b * b);
    }

    float D(const Vector3f& wh) const {
        float tan2Theta = std::pow(Frame::tanTheta(wh), 2);
        float cos3Theta = std::pow(Frame::cosTheta(wh), 3);
        if (std::isinf(tan2Theta) || cos3Theta <= 0) {
            return 0.0f;
        }
        return INV_PI * std::exp(-tan2Theta / m_alpha2) / (m_alpha2 * cos3Theta);
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord& bRec) const {
        if (bRec.measure != ESolidAngle) {
            return Color3f(0.0f);
        }

        // fresnel term
        float F = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
        float cosThetaI = Frame::cosTheta(bRec.wi);
        float cosThetaO = Frame::cosTheta(bRec.wo);
        if (cosThetaI == 0 || cosThetaO == 0) {
            return Color3f(0.0f);
        }

        Color3f diffuseTerm(0.0f), reflectanceTerm(0.0f), transmissionTerm(0.0f);

        bool sameSide = bRec.wi.z() * bRec.wo.z() > 0;

        if (sameSide) {
            Vector3f wh = (bRec.wi + bRec.wo).normalized() * (bRec.wi.z() < 0 ? -1 : 1);
            float cosThetaH = Frame::cosTheta(wh);
            // distribution term
            float D = this->D(wh);
            // geometry term
            float G = this->G1(bRec.wi, wh) * this->G1(bRec.wo, wh);
            reflectanceTerm = m_ks * D * F * G / (4 * cosThetaI * cosThetaH * cosThetaO);
            diffuseTerm = m_kd * INV_PI;
        }
        else {
            float eta = cosThetaI > 0 ? (m_intIOR / m_extIOR) : (m_extIOR / m_intIOR);
            Vector3f wh = (bRec.wi + bRec.wo * eta).normalized();
            if (wh.z() < 0) {
                wh = -wh;
            }
            float cosThetaH = Frame::cosTheta(wh);
            float sqrtDenom = wh.dot(bRec.wi) + wh.dot(bRec.wo) * eta;
            // distribution term
            float D = this->D(wh);
            // geometry term
            float G = this->G1(bRec.wi, wh) * this->G1(bRec.wo, wh);
            transmissionTerm = m_ks * (1.0f - F) * std::abs(D * G * eta * eta * std::abs(wh.dot(bRec.wi)) * std::abs(wh.dot(bRec.wo)) /
                (sqrtDenom * sqrtDenom * cosThetaI * cosThetaO * cosThetaH));
        }

        Color3f L = diffuseTerm + reflectanceTerm + transmissionTerm;
        return L;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord& bRec) const {
        float cosThetaO = Frame::cosTheta(bRec.wo);
        float cosThetaI = Frame::cosTheta(bRec.wi);

        bool sameSide = bRec.wi.z() * bRec.wo.z() > 0;
        float pdf = 0.0f;
        if (sameSide) {
            Vector3f wh = (bRec.wi + bRec.wo).normalized() * (bRec.wi.z() < 0 ? -1 : 1);

            float jacobian = 1.0f / (4 * std::abs(wh.dot(bRec.wo)));
            pdf = m_ks * this->D(wh) * jacobian + (1 - m_ks) * INV_PI * cosThetaO;
        }
        else {
            float eta = cosThetaI > 0 ? (m_intIOR / m_extIOR) : (m_extIOR / m_intIOR);
            Vector3f wh = (bRec.wi + bRec.wo * eta).normalized();
            if (wh.z() < 0) {
                wh = -wh;
            }
            float sqrtDenom = wh.dot(bRec.wi) + wh.dot(bRec.wo) * eta;
            float jacobian = std::abs(eta * eta * wh.dot(bRec.wo) / (sqrtDenom * sqrtDenom));
            pdf = m_ks * this->D(wh) * jacobian;
        }

        return pdf;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord& bRec, const Point2f& _sample) const {
        Point2f sample(_sample);
        bRec.measure = ESolidAngle;

        if (sample.x() < m_ks) {
            // reuse sample.x()
            sample.x() /= m_ks;
            
            float F = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
            if (sample.x() < F) {
                sample.x() /= F;
                Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
                bRec.wo = reflect(bRec.wi, wh);
                bRec.wh = wh;
            }
            else {
                sample.x() = (sample.x() - F) / (1.0f - F);
                Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
                bool entering = Frame::cosTheta(bRec.wi) > 0;
                bRec.eta = entering ? (m_intIOR / m_extIOR) : (m_extIOR / m_intIOR);
                bRec.wh = wh;
                if (!refract(bRec.wi, wh, 1.0f / bRec.eta, bRec.wo)) {
                    return Color3f(0.0f);
                }
            }
        }
        else {
            sample.x() = (sample.x() - m_ks) / (1.0f - m_ks);
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }

        float pdf = this->pdf(bRec);
        if (pdf <= 0) {
            return Color3f(0.0f);
        }
        return this->eval(bRec) * std::abs(Frame::cosTheta(bRec.wo)) / pdf;

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "RoughDieletric[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha, m_alpha2;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(RoughDieletric, "rough_dielectric");
NORI_NAMESPACE_END
