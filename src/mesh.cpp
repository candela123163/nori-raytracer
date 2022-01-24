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

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>
#include <nori/accel.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

Color3f Intersection::Le(const Vector3f& wi) const {
    const Emitter* area = mesh->getEmitter();
    return area ? area->Li(*this, wi) : Color3f(0.0f);
}

bool VisibilityTester::Unoccluded(const Scene& scene) const {
    return !scene.rayIntersect(p0.SpawnRayTo(p1));
}

Mesh::Mesh() { 
    m_accel = std::make_unique<Accel>();
}

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }
    m_accel->addMesh(this);
    m_accel->build();
    this->initDPdf();
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

Intersection Mesh::sample(const Intersection& ref, const Point2f& sample, float& pdf) const {
    Intersection its;
    Point2f sample_(sample);
    // sample a triangle
    float trianglePdf = 0.0f;
    size_t triangleIdx = m_dpdf.sampleReuse(sample_.x(), trianglePdf);
    
    // sample on area
    float sqrtTerm = std::sqrt(1 - sample_.x());
    float alpha = 1 - sqrtTerm;
    float beta = sample_.y() * sqrtTerm;
    Vector3f bary(1 - alpha - beta, alpha, beta);

    // build intersection data
    uint32_t idx0 = m_F(0, triangleIdx), idx1 = m_F(1, triangleIdx), idx2 = m_F(2, triangleIdx);
    Point3f p0 = m_V.col(idx0), p1 = m_V.col(idx1), p2 = m_V.col(idx2);
    its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;
    its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());
    if (m_N.size() > 0) {
        its.shFrame = Frame(
            (bary.x() * m_N.col(idx0) +
             bary.y() * m_N.col(idx1) +
             bary.z() * m_N.col(idx2)).normalized());
    }
    else {
        its.shFrame = its.geoFrame;
    }
    its.mesh = this;
    its.triangleIdx = triangleIdx;
    
    // calculate pdf
    float areaPdf = 1.0f / this->surfaceArea(triangleIdx);
    pdf = trianglePdf * areaPdf;
    // convert from area measure to solid angle measure
    Vector3f wi = its.p - ref.p;
    float squaredDistance = wi.squaredNorm();
    wi.normalize();
    its.wi = its.toLocal(-wi);
    pdf *= squaredDistance / its.shFrame.cosTheta(its.wi);
    if (std::isinf(pdf) || pdf <= 0.0f) {
        pdf = 0.0f;
    }

    return its;
}

float Mesh::pdf(const Intersection& ref, const Vector3f& wo) const {
    Ray3f ray = ref.SpawnRay(wo);
    Intersection its;
    if (!m_accel->rayIntersect(ray, its, false)) {
        return 0.0f;
    }
    float trianglePdf = m_dpdf[its.triangleIdx];
    float areaPdf = 1.0f / this->surfaceArea(its.triangleIdx);
    float pdf = trianglePdf * areaPdf;
    Vector3f localWi = its.shFrame.toLocal(-wo);
    pdf *= (ref.p - its.p).squaredNorm() / its.shFrame.cosTheta(localWi);
    if (std::isinf(pdf)  || pdf <= 0.0f) {
        pdf = 0.0f;
    }
    return pdf;
}

void Mesh::initDPdf() {
    m_dpdf.clear();
    m_dpdf.reserve(this->getTriangleCount());
    for (size_t i = 0; i < this->getTriangleCount(); i++)
    {
        m_dpdf.append(this->surfaceArea(i));
    }
    m_dpdf.normalize();
}

NORI_NAMESPACE_END
