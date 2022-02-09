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

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
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
        if (Frame::cosTheta(wv) <= 0 || wv.dot(wh) <= 0) {
            return 0.0f;
        }
        float b = 1.0f / (m_alpha * Frame::tanTheta(wv));
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
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle) {
            return Color3f(0.0f);
        }

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float cosThetaI = Frame::cosTheta(bRec.wi);
        float cosThetaH = Frame::cosTheta(wh);
        float cosThetaO = Frame::cosTheta(bRec.wo);
        if (cosThetaI <= 0 || cosThetaH <= 0 || cosThetaO <= 0) {
            return Color3f(0.0f);
        }

        // distribution term
        float D = this->D(wh);
        
        // fresnel term
        float F = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

        // geometry term
        float G = this->G1(bRec.wi, wh) * this->G1(bRec.wo, wh);
        
        // diffuse term
        Color3f diffuse_term = m_kd * INV_PI;

        Color3f L = diffuse_term + m_ks * D * F * G / (4 * cosThetaI * cosThetaH * cosThetaO);

        return L;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        float cosThetaO = Frame::cosTheta(bRec.wo);
        if (cosThetaO <= 0) {
            return 0.0f;
        }
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float jacobian = 1.0f / (4 * wh.dot(bRec.wo));
        return m_ks * this->D(wh) * jacobian + (1 - m_ks) * INV_PI * cosThetaO;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        Point2f sample(_sample);
        
        if (sample.x() < m_ks) {
            // reuse sample.x()
            sample.x() /= m_ks;
            Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
            bRec.wo = reflect(bRec.wi, wh);
        }
        else {
            sample.x() = (sample.x() - m_ks) / (1.0f - m_ks);
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }

        bRec.measure = ESolidAngle;        
        float pdf = this->pdf(bRec);
        if (pdf <= 0) {
            return Color3f(0.0f);
        }
        return this->eval(bRec) * std::max(0.0f, Frame::cosTheta(bRec.wo)) / pdf;

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
            "Microfacet[\n"
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

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
