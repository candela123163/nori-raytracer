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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    float u[2] = { sample.x(), sample.y() };
    float x[2];
    for (size_t i = 0; i < 2; i++)
    {
        if (u[i] <= 0.5) {
            x[i] = -1 + std::sqrt(2 * u[i]);
        }
        else {
            x[i] = 1 - std::sqrt(2 - 2 * u[i]);
        }
    }
    return { x[0], x[1] };
}

float Warp::squareToTentPdf(const Point2f &p) {
    return (1 - std::abs(p.x())) * (1 - std::abs(p.y()));
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = std::sqrt(sample.x());
    float theta = 2 * M_PI * sample.y();
    return { r * std::cos(theta), r * std::sin(theta) };
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return p.squaredNorm() <= 1.0f ? INV_PI : 0.0f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float cosTheta = 1 - 2 * sample.x();
    float sinTheta = std::sqrt(1 - cosTheta * cosTheta);
    float phi = 2 * M_PI * sample.y();
    return { std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta };
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return INV_FOURPI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float cosTheta = sample.x();
    float sinTheta = std::sqrt(1 - cosTheta * cosTheta);
    float phi = 2 * M_PI * sample.y();
    return { std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta };
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return v.z() >= 0.0f ? INV_TWOPI : 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Point2f d = Warp::squareToUniformDisk(sample);
    float z = std::sqrt(1 - d.x() * d.x() - d.y() * d.y());
    return { d.x(), d.y(), z };
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    float cosTheta = v.z();
    return v.z() >= 0.0f ? cosTheta * INV_PI : 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float phi = 2 * M_PI * sample.x();
    float tan2Theta = -alpha * alpha * std::log(sample.y());
    float cosTheta = 1.0f / std::sqrt(1.0f + tan2Theta);
    float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);
    return { sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta };
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    float cosTheta = m.z();
    float cos3Theta = cosTheta * cosTheta * cosTheta;
    float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);
    float tan2Theta = std::pow(sinTheta / cosTheta, 2);
    float alpha2 = alpha * alpha;
    float pdf = INV_PI * std::exp(-tan2Theta / alpha2) / (alpha2 * cos3Theta);
    return m.z() > 0.0f ? pdf : 0.0f;
}

NORI_NAMESPACE_END
