#ifndef SIMPLEX_MESH_H
#define SIMPLEX_MESH_H

#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/DataStructure/HashTable.h>
#include <vector>

namespace ZIRAN {
/**
    This is a class for keeping track of simplex mesh topology: e, i, j, k (element e has nodes i,j,k), etc. The class also has basic helpers like routines for computing the boundary simplex mesh
*/
template <int dim>
class SimplexMesh {
public:
    static constexpr int manifold_dim = dim;
    static constexpr int num_vertices = dim + 1;
    static constexpr int num_faces = dim + 1;
    static constexpr int num_edges = 3 * dim - 3;
    static constexpr int boundary_vertices = dim;
    typedef Eigen::Vector2i NIV;
    typedef Vector<int, boundary_vertices> BIV;
    typedef Vector<int, 2> EIV;
    typedef Vector<int, num_vertices> IV;
    typedef Vector<int, num_edges> EEIV; // Element Edge Indices Vector
    typedef std::array<BIV, num_faces> Faces;
    typedef Eigen::SparseMatrix<int, Eigen::RowMajor> SpMat;
    StdVector<IV> indices;
    StdVector<BIV> boundary_indices;
    SpMat connectivity;
    bool neighboring_initialized = false;

    // Mesh Class owned neighboring infos
    bool is_neighbor_list_created = false;
    bool neighbor_list_sorted;
    StdVector<NIV> neighbor_list;
    StdVector<BIV> neighbor_vertex_list;

    HashTable<int, StdVector<int>> vertex_to_neighbor_vertices;
    bool vertex_to_neighbor_vertices_initialized = false;

    HashTable<int, StdVector<int>> vertex_to_neighbor_elements;
    bool vertex_to_neighbor_elements_initialized = false;

    SimplexMesh()
    {
    }

    SimplexMesh(StdVector<IV>&& indices)
        : indices(indices)
    {
    }

    inline int index(const size_t e, const int i) const { return indices[e](i); }

    inline size_t numberElements() const { return indices.size(); }

    /**
        Find boundary elements for the mesh using hash table
        Any (dim -1) dimension element (which is stored in vector of size dim) traversed only once is a boundary element.
        Any (dim -1) dimension element traversed twice is an interior element.
    */

    void initializeBoundaryElements();

    void constructListOfEdges(StdVector<EIV>& edge_list, HashTable<EIV, int>* edge_hash_ptr = nullptr);

    // void constructListOfEdges(StdVector<EIV>& edge_list, EEIV* associated_edge_list_begin, EEIV* associated_edge_list_end); // TODO: this is no longer needed I believe (I == Stephanie)

    void constructListOfFaces(StdVector<BIV>& face_list, HashTable<BIV, int>* face_hash_ptr = nullptr);

    void getFaces(const size_t element, std::array<Vector<int, dim>, dim + 1>& result);

    //return the two/three vertices of the common face, sorted w.r.t. vertex ID
    static void findSharedFace(const IV& element1, const IV& element2, BIV& face);

    static void findSharedFaceAndIndex(const IV& element1, const IV& element2, BIV& sorted, BIV* local_index1 = nullptr, BIV* local_index2 = nullptr);

    static void findUnsortedFace(const IV& element, const BIV& sorted, BIV& face);

    // construct list of pair of endpoints for the interior edges
    void constructNeighborVertexList(const StdVector<NIV>& neighbor_list, StdVector<BIV>& neighbor_vertex_list) const;

    void constructNeighborVertexListInternal();

    static const char* name();

    /**
        Construct a list of face touching elements.
        Optional sorting of each element pair.
    */
    void constructListOfFaceNeighboringElements(StdVector<Vector<int, 2>>& result, const bool sort = false, HashTable<BIV, Eigen::Vector2i>* face_hash_ptr = nullptr);

    void constructListOfFaceNeighboringElementsInternal(const bool sort = false, HashTable<BIV, Eigen::Vector2i>* face_hash_ptr = nullptr);

    void initializeNeighboringElements();

    // self is not considered as a neighbor
    void initializeVertexToNeighbouringVertices();

    void initializeVertexToNeighbouringElements();

    /**
       Input: index is an element index.
       Input: target is a list of element indicies in the same mesh.
       Output: result is the closest one in the target to the index.
               where the distance is judged by graph connectivity.
    */
    void closestElement(const int index, int& result, const HashTable<int, bool>& target);

    void findIndexBounds(int& min, int& max) const;
};
} // namespace ZIRAN

#endif
