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

#pragma once

#include <nori/mesh.h>
#include <memory>

NORI_NAMESPACE_BEGIN

/*
* BVH data structure declaratins
*/
struct BVHBuildNode
{
    std::unique_ptr<BVHBuildNode> left, right;
    BoundingBox3f bbox;
    uint32_t firstPrimOffset = 0, nPrimitives = 0;

    void setLeaf(uint32_t firstPrimOffset_, uint32_t nPrimitives_, BoundingBox3f bbox_) {
        left = nullptr;
        right = nullptr;
        firstPrimOffset = firstPrimOffset_;
        nPrimitives = nPrimitives_;
        bbox = bbox_;
    }

    void setInterior(std::unique_ptr<BVHBuildNode>&& left_, std::unique_ptr<BVHBuildNode>&& right_) {
        left = std::move(left_);
        right = std::move(right_);
        bbox.expandBy(left->bbox);
        bbox.expandBy(right->bbox);
    }
};

struct FlatBVHNode
{
    BoundingBox3f bbox;
    union {
        uint32_t primitivesOffset;
        uint32_t secondChildOffset;
    };
    uint32_t nPrimitives;   // 0 -> interior node
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
    Accel();
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    /*
    * \brief build bvh tree recursively
    * return root node of bvh tree
    * 
    * \param primitiveIndices
    *     Indices of primiitives
    * 
    * \param start
    *     First index of primitive of this node
    * 
    * \param end 
    *     1 + Last index of primitive of this node
    * 
    * \param nodeCount
    *     the number of nodes
    */
    std::unique_ptr<BVHBuildNode> recursiveBuild(std::vector<uint32_t>& primitiveIndices, uint32_t start, uint32_t end, uint32_t& nodeCount);

    /*
    * \brief convert bvh binary tree to linear array for faster search
    * return index of primitive of this node in ordered indices
    * 
    * \param node
    *     Current BVHBuildNode
    * 
    * \param nodeIdx
    *     Current BVHBuildNode's index
    */
    uint32_t flattenBVHTree(const std::unique_ptr<BVHBuildNode>& node, uint32_t& nodeIdx);

    /*
    * \brief query mesh by global triangle index
    * return pointer to result mesh
    *
    *
    * \param globalTriangleIdx
    *     global triangle index
    * 
    * \param localTriangleIdx
    *     triangle index of queried mesh
    */
    Mesh* queryMesh(size_t& last_left_bound, size_t& last_right_bound, size_t globalTriangleIdx, size_t& localTriangleIdx) const;

private:
    //Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    std::vector<Mesh*> m_meshes;
    std::vector<size_t> m_mesh_triangles;   ///< cumulative mesh triangle number(i.e. m_mesh_triangles[0] = 0, m_mesh_triangles[1] = first mesh triangle number
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
    std::vector<FlatBVHNode> m_BVH_nodes; ///< array of bvh nodes
    std::vector<uint32_t> m_ordered_indices; ///< ordered indices of mesh primitives
};

NORI_NAMESPACE_END
