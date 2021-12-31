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

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <stack>

NORI_NAMESPACE_BEGIN

#define BVH_LEAF_MAX_PRIMITIVES 10


std::unique_ptr<BVHBuildNode> Accel::recursiveBuild(std::vector<uint32_t>& primitiveIndices, uint32_t start, uint32_t end, uint32_t& nodeCount) {
    ++nodeCount;
    auto build_node = std::make_unique<BVHBuildNode>();
    
    BoundingBox3f bbox;
    for (uint32_t i = start; i < end; i++)
    {
        uint32_t index = primitiveIndices[i];
        bbox.expandBy(m_mesh->getBoundingBox(index));
    }

    int nPrimitives = end - start;
    // leaf node
    if (nPrimitives <= BVH_LEAF_MAX_PRIMITIVES || !bbox.isValid()) {
        int firstPrimOffset = m_ordered_indices.size();
        for (uint32_t i = start; i < end; i++) {
            uint32_t index = primitiveIndices[i];
            m_ordered_indices.push_back(index);
        }
        build_node->setLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    // interior node
    else {
        // try partition with mid pos of boundingbox's largest axis
        int largestAxis = bbox.getLargestAxis();
        float axisMidPos = bbox.getCenter()(largestAxis);
        auto iter = std::partition(
            primitiveIndices.begin() + start,
            primitiveIndices.begin() + end,
            [axisMidPos, largestAxis, this](uint32_t _idx) {
            return this->m_mesh->getBoundingBox(_idx).getCenter()(largestAxis) < axisMidPos;
        });
        int mid = iter - primitiveIndices.begin();

        // partition failed, partition primitives into equally-sized subsets
        if (mid == start || mid == end) {
            mid = (start + end) / 2;
            std::nth_element(
                primitiveIndices.begin() + start,
                primitiveIndices.begin() + mid,
                primitiveIndices.begin() + end,
                [largestAxis, this](uint32_t _left, uint32_t _right) {
                return this->m_mesh->getBoundingBox(_left).getCenter()(largestAxis) < this->m_mesh->getBoundingBox(_right).getCenter()(largestAxis);
            });
        }

        build_node->setInterior(this->recursiveBuild(primitiveIndices, start, mid, nodeCount),
            this->recursiveBuild(primitiveIndices, mid, end, nodeCount));
    }
    return build_node;
}

uint32_t Accel::flattenBVHTree(const std::unique_ptr<BVHBuildNode>& node, uint32_t& nodeIdx) {
    auto& currentNode = m_BVH_nodes[nodeIdx];
    currentNode.bbox = node->bbox;
    int currentIdx = nodeIdx++;
    // leaf node
    if (node->nPrimitives > 0) {
        currentNode.primitivesOffset = node->firstPrimOffset;
        currentNode.nPrimitives = node->nPrimitives;
    }
    // interior node
    else {
        this->flattenBVHTree(node->left, nodeIdx);
        currentNode.secondChildOffset = this->flattenBVHTree(node->right, nodeIdx);
    }
    return currentIdx;
}

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
    m_ordered_indices.clear();
    std::vector<uint32_t> primitiveIndices(m_mesh->getTriangleCount(), 0);
    for (uint32_t i = 0; i < primitiveIndices.size(); i++)
    {
        primitiveIndices[i] = i;
    }
    uint32_t bvhNodeCount = 0;
    uint32_t flatBVHNodeIdx = 0;
    auto root = this->recursiveBuild(primitiveIndices, 0, primitiveIndices.size(), bvhNodeCount);
    m_BVH_nodes.resize(bvhNodeCount);
    this->flattenBVHTree(root, flatBVHNodeIdx);
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Search in flatten bvh tree */
    std::stack<uint32_t> travelStack;

    /* push root index into stack */
    travelStack.push(0);
    while (!travelStack.empty()) {
        uint32_t curIdx = travelStack.top();
        travelStack.pop();
        const auto& curNode = m_BVH_nodes[curIdx];
        
        if (!curNode.bbox.rayIntersect(ray)) {
            continue;
        }

        if (curNode.nPrimitives > 0) {
            for (uint32_t i = 0; i < curNode.nPrimitives; i++) {
                uint32_t primitiveIdx = m_ordered_indices[curNode.primitivesOffset + i];
                float u, v, t;
                if (m_mesh->rayIntersect(primitiveIdx, ray, u, v, t)) {
                    /* An intersection was found! Can terminate
                       immediately if this is a shadow ray query */
                    if (shadowRay)
                        return true;
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    f = primitiveIdx;
                    foundIntersection = true;
                }
            }
        }
        else {
            uint32_t leftChildIdx = curIdx + 1;
            uint32_t rightChildIdx = curNode.secondChildOffset;

            float leftNearT, leftFarT, rightNearT, rightFarT;
            bool leftIntersect = m_BVH_nodes[leftChildIdx].bbox.rayIntersect(ray, leftNearT, leftFarT);
            bool rightIntersect = m_BVH_nodes[rightChildIdx].bbox.rayIntersect(ray, rightNearT, rightFarT);

            if (leftIntersect && rightIntersect) {
                if (leftNearT < rightNearT)
                {
                    travelStack.push(rightChildIdx);
                    travelStack.push(leftChildIdx);
                }
                else {
                    travelStack.push(leftChildIdx);
                    travelStack.push(rightChildIdx);
                }
            }
            else {
                if (leftIntersect) {
                    travelStack.push(leftChildIdx);
                }
                if (rightIntersect) {
                    travelStack.push(rightChildIdx);
                }
            }
        }
    }


    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

