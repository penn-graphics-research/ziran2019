#include <Ziran/Math/Geometry/SimplexMesh.h>
#include <Ziran/CS/DataStructure/HashTable.h>
#include <Ziran/CS/Util/Debug.h>
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <queue>

namespace ZIRAN {

template <int dim>
void SimplexMesh<dim>::initializeBoundaryElements()
{
    HashTable<BIV, BIV> face_hash;
    size_t number_elements = indices.size();
    Faces faces;

    for (size_t i = 0; i < number_elements; i++) {
        //go through each element
        //and get all faces for each element
        getFaces(i, faces);

        //for each face, sort indices
        //if the set of indices does not exist, put in hash table
        // if the indices already exists, delete from hash table
        // all faces left are boundary faces
        for (size_t j = 0; j < faces.size(); j++) {
            BIV sorted = faces[j];
            std::sort(sorted.data(), sorted.data() + sorted.size());

            if (!face_hash.erase(sorted))
                face_hash.insert(sorted, faces[j]);
        }
    }
    boundary_indices.clear();
    for (auto it = face_hash.cbegin(); it != face_hash.cend(); it++) {
        boundary_indices.emplace_back(std::move(it->value));
    }
}

template <>
void SimplexMesh<1>::getFaces(
    const size_t element,
    std::array<Vector<int, 1>, 2>& result)
{
    int n0 = indices[element](0);
    int n1 = indices[element](1);
    result[0](0) = n0;
    result[1](0) = n1;
}

template <>
void SimplexMesh<2>::getFaces(
    const size_t element,
    std::array<Vector<int, 2>, 3>& result)
{
    int n0 = indices[element](0);
    int n1 = indices[element](1);
    int n2 = indices[element](2);
    result[0] = Vector<int, 2>(n0, n1);
    result[1] = Vector<int, 2>(n1, n2);
    result[2] = Vector<int, 2>(n2, n0);
}

template <>
void SimplexMesh<3>::getFaces(
    const size_t element,
    std::array<Vector<int, 3>, 4>& result)
{
    int n0 = indices[element](0);
    int n1 = indices[element](1);
    int n2 = indices[element](2);
    int n3 = indices[element](3);
    result[0] = Vector<int, 3>(n0, n1, n2);
    result[1] = Vector<int, 3>(n1, n2, n3);
    result[2] = Vector<int, 3>(n2, n3, n0);
    result[3] = Vector<int, 3>(n3, n0, n1);
}

template <>
void SimplexMesh<1>::constructListOfEdges(StdVector<EIV>& edge_list, HashTable<EIV, int>* edge_hash_ptr)
{

    HashTable<EIV, int> inner_edge_hash;
    if (!edge_hash_ptr)
        edge_hash_ptr = &inner_edge_hash;

    size_t number_elements = indices.size();
    Faces faces;

    edge_list.clear();

    EIV sorted;
    for (size_t i = 0; i < number_elements; ++i) {
        sorted = indices[i];
        std::sort(sorted.data(), sorted.data() + sorted.size());
        auto result = edge_hash_ptr->insert(sorted, (int)edge_list.size());

        if (result.inserted)
            edge_list.emplace_back(sorted);
    }
}

template <>
void SimplexMesh<2>::constructListOfEdges(StdVector<EIV>& edge_list, HashTable<EIV, int>* edge_hash_ptr)
{
    SimplexMesh<2>::constructListOfFaces(edge_list, edge_hash_ptr);
}

template <>
void SimplexMesh<3>::constructListOfEdges(StdVector<EIV>& edge_list, HashTable<EIV, int>* edge_hash_ptr)
{
    HashTable<EIV, int> inner_edge_hash;
    if (!edge_hash_ptr)
        edge_hash_ptr = &inner_edge_hash;

    size_t number_elements = indices.size();
    Faces faces;

    edge_list.clear();

    for (size_t i = 0; i < number_elements; i++) {
        getFaces(i, faces);

        for (size_t j = 0; j < faces.size(); j++) {
            EIV edge;
            for (int k = 0; k < 3; k++) {
                edge(0) = faces[j](k);
                edge(1) = faces[j]((k + 1) % 3);
                EIV sorted = edge;
                std::sort(sorted.data(), sorted.data() + sorted.size());
                auto result = edge_hash_ptr->insert(sorted, (int)edge_list.size());
                if (result.inserted)
                    edge_list.emplace_back(sorted);
            }
        }
    }
}

// template <>
// void SimplexMesh<1>::constructListOfEdges(StdVector<EIV>& edge_list, EEIV* associated_edge_list_begin, EEIV* associated_edge_list_end)
// {
//     // associated edge list doesn't make sense for seg mesh
//     constructListOfEdges(edge_list);
// }

// template <>
// void SimplexMesh<2>::constructListOfEdges(StdVector<EIV>& edge_list, EEIV* associated_edge_list_begin, EEIV* associated_edge_list_end)
// {
//     HashTable<BIV, size_t> face_hash;
//     size_t number_elements = indices.size();
//     Faces faces;

//     for (size_t i = 0; i < number_elements; i++) {
//         getFaces(i, faces);
//         for (size_t j = 0; j < faces.size(); j++) {
//             BIV sorted = faces[j];
//             std::sort(sorted.data(), sorted.data() + sorted.size());
//             auto result = face_hash.insert(sorted, edge_list.size());

//             if (result.inserted)
//                 edge_list.emplace_back(sorted);
//             assert(associated_edge_list_begin + i < associated_edge_list_end);
//             associated_edge_list_begin[i](j) = result.value;
//         }
//     }
// }

// template <>
// void SimplexMesh<3>::constructListOfEdges(StdVector<EIV>& edge_list, EEIV* associated_edge_list_begin, EEIV* associated_edge_list_end)
// {
//     //TODO
//     constructListOfEdges(edge_list);
// }

template <int dim>
void SimplexMesh<dim>::constructListOfFaces(StdVector<BIV>& face_list, HashTable<BIV, int>* face_hash_ptr)
{
    HashTable<BIV, int> inner_face_hash;
    if (!face_hash_ptr)
        face_hash_ptr = &inner_face_hash;

    size_t number_elements = indices.size();
    Faces faces;

    face_list.clear();

    for (size_t i = 0; i < number_elements; i++) {
        getFaces(i, faces);
        for (size_t j = 0; j < faces.size(); j++) {
            BIV sorted = faces[j];
            std::sort(sorted.data(), sorted.data() + sorted.size());
            auto result = face_hash_ptr->insert(sorted, (int)face_list.size());
            if (result.inserted)
                face_list.emplace_back(sorted);
        }
    }
}

template <int dim>
void SimplexMesh<dim>::findSharedFace(const IV& element1, const IV& element2, BIV& face)
{
    int e = 0;
    bool job_done = false;
    for (int i = 0; i < num_vertices && !job_done; i++) {
        for (int j = 0; j < num_vertices && !job_done; j++) {
            if (e == boundary_vertices)
                job_done = true;
            if (element1(i) == element2(j) && !job_done) {
                face(e) = element1(i);
                e++;
            }
        }
    }

    std::sort(face.data(), face.data() + face.size());
}

template <int dim>
void SimplexMesh<dim>::findSharedFaceAndIndex(const IV& element1, const IV& element2, BIV& sorted, BIV* local_index1_ptr, BIV* local_index2_ptr)
{
    BIV inner_local_index1, inner_local_index2;
    if (!local_index1_ptr)
        local_index1_ptr = &inner_local_index1;
    if (!local_index2_ptr)
        local_index2_ptr = &inner_local_index2;

    int find = 0;
    for (int k = 0; k < num_vertices && find < boundary_vertices; k++) {
        for (int l = 0; l < num_vertices && find < boundary_vertices; l++) {
            if (element1(k) == element2(l)) {
                (*local_index1_ptr)(find) = k;
                (*local_index2_ptr)(find) = l;
                sorted(find) = element1(k);
                find++;
                break;
            }
        }
    }
    std::sort(local_index1_ptr->data(), local_index1_ptr->data() + local_index1_ptr->size());
    std::sort(local_index2_ptr->data(), local_index2_ptr->data() + local_index2_ptr->size());
    std::sort(sorted.data(), sorted.data() + sorted.size());
}

template <int dim>
void SimplexMesh<dim>::findUnsortedFace(const IV& element, const BIV& sorted, BIV& face)
{
    int find = 0;
    for (int k = 0; k < num_vertices && find < boundary_vertices; k++) {
        for (int l = 0; l < boundary_vertices && find < boundary_vertices; l++) {
            if (element(k) == sorted(l)) {
                face(find) = element(k);
                find++;
                break;
            }
        }
    }
}

template <int dim>
void SimplexMesh<dim>::constructNeighborVertexList(const StdVector<NIV>& neighbor_list, StdVector<BIV>& neighbor_vertex_list) const
{
    neighbor_vertex_list.clear();
    neighbor_vertex_list.reserve(neighbor_list.size());
    for (int n = 0; n < int(neighbor_list.size()); n++) {
        NIV face = neighbor_list[n];
        IV element1 = indices[face(0)];
        IV element2 = indices[face(1)];
        BIV vertices;
        findSharedFace(element1, element2, vertices);
        neighbor_vertex_list.emplace_back(vertices);
    }
}

template <int dim>
void SimplexMesh<dim>::constructNeighborVertexListInternal()
{
    // if the neighbor_list hasn't been created, then create one.
    if (!is_neighbor_list_created) {
        constructListOfFaceNeighboringElementsInternal(false);
    }

    neighbor_vertex_list.clear();
    neighbor_vertex_list.reserve(neighbor_list.size());
    for (int n = 0; n < int(neighbor_list.size()); n++) {
        NIV face = neighbor_list[n];
        IV element1 = indices[face(0)];
        IV element2 = indices[face(1)];
        BIV vertices;
        findSharedFace(element1, element2, vertices);
        neighbor_vertex_list.emplace_back(vertices);
    }
}

template <int dim>
const char* SimplexMesh<dim>::name()
{
    switch (dim) {
    case 1:
        return "Seg";
        break;
    case 2:
        return "Tri";
        break;
    case 3:
        return "Tet";
        break;
    default:
        ZIRAN_ASSERT("Unsupported Element dimension ", dim);
    }
}

/**
        Construct a list of face touching elements.
        Optional sorting of each element pair.
    */
template <int dim>
void SimplexMesh<dim>::constructListOfFaceNeighboringElements(StdVector<Eigen::Vector2i>& result, const bool sort, HashTable<BIV, Eigen::Vector2i>* face_hash_ptr)
{
    HashTable<BIV, Eigen::Vector2i> inner_face_hash;
    if (!face_hash_ptr)
        face_hash_ptr = &inner_face_hash;

    size_t number_elements = indices.size();
    Faces faces;

    for (size_t i = 0; i < number_elements; i++) {
        getFaces(i, faces);
        for (size_t j = 0; j < faces.size(); j++) {
            BIV sorted = faces[j];
            std::sort(sorted.data(), sorted.data() + sorted.size());

            Eigen::Vector2i* v = face_hash_ptr->get(sorted);
            if (v == nullptr)
                face_hash_ptr->insert(sorted, Eigen::Vector2i(i, -1));
            else {
                ZIRAN_ASSERT((*v)(0) != -1 && (*v)(0) != (int)i && (*v)(1) == -1);
                (*v)(1) = (int)i;
            }
        }
    }

    result.clear();
    for (auto it = face_hash_ptr->cbegin(); it != face_hash_ptr->cend(); it++) {
        Eigen::Vector2i elements = it->value.template head<2>();
        ZIRAN_ASSERT(elements(0) != -1);
        if (elements(1) != -1)
            result.push_back(elements);
    }

    if (sort) {
        for (size_t i = 0; i < result.size(); i++) {
            std::sort(result[i].data(), result[i].data() + result[i].size());
        }
    }
}

template <int dim>
void SimplexMesh<dim>::constructListOfFaceNeighboringElementsInternal(const bool sort, HashTable<BIV, Eigen::Vector2i>* face_hash_ptr)
{
    neighbor_list_sorted = sort;
    HashTable<BIV, Eigen::Vector2i> inner_face_hash;
    if (!face_hash_ptr)
        face_hash_ptr = &inner_face_hash;

    size_t number_elements = indices.size();
    Faces faces;

    for (size_t i = 0; i < number_elements; i++) {
        getFaces(i, faces);
        for (size_t j = 0; j < faces.size(); j++) {
            BIV sorted = faces[j];
            std::sort(sorted.data(), sorted.data() + sorted.size());

            Eigen::Vector2i* v = face_hash_ptr->get(sorted);
            if (v == nullptr)
                face_hash_ptr->insert(sorted, Eigen::Vector2i(i, -1));
            else {
                ZIRAN_ASSERT((*v)(0) != -1 && (*v)(0) != (int)i && (*v)(1) == -1);
                (*v)(1) = (int)i;
            }
        }
    }

    neighbor_list.clear();
    for (auto it = face_hash_ptr->cbegin(); it != face_hash_ptr->cend(); it++) {
        Eigen::Vector2i elements = it->value.template head<2>();
        ZIRAN_ASSERT(elements(0) != -1);
        if (elements(1) != -1)
            neighbor_list.push_back(elements);
    }

    if (sort) {
        for (size_t i = 0; i < neighbor_list.size(); i++) {
            std::sort(neighbor_list[i].data(), neighbor_list[i].data() + neighbor_list[i].size());
        }
    }

    // set flag recording that neighbor_list has been created
    is_neighbor_list_created = true;
}

template <int dim>
void SimplexMesh<dim>::initializeNeighboringElements()
{
    StdVector<Eigen::Vector2i> result;
    constructListOfFaceNeighboringElements(result, false);
    const int num_elements = indices.size();
    connectivity.resize(num_elements, num_elements);
    for (size_t i = 0; i < result.size(); i++) {
        Eigen::Vector2i pair = result[i];
        connectivity.insert(pair(0), pair(1)) = 1;
        connectivity.insert(pair(1), pair(0)) = 1;
    }
    neighboring_initialized = true;
}

// self is not considered as a neighbor
template <int dim>
void SimplexMesh<dim>::initializeVertexToNeighbouringVertices()
{
    if (!vertex_to_neighbor_vertices_initialized) {
        size_t number_elements = indices.size();
        Faces faces;
        for (size_t i = 0; i < number_elements; i++) {
            for (int p = 0; p < dim + 1; p++) {
                int my_index = indices[i](p);
                for (int q = 0; q < dim + 1; q++) {
                    int your_index = indices[i](q);
                    if (p != q) {
                        StdVector<int>& list = vertex_to_neighbor_vertices[my_index];
                        auto found = std::find(std::begin(list), std::end(list), your_index);
                        if (found == list.end())
                            list.emplace_back(your_index);
                    }
                }
            }
        }
        vertex_to_neighbor_vertices_initialized = true;
    }
}

template <int dim>
void SimplexMesh<dim>::initializeVertexToNeighbouringElements()
{
    if (!vertex_to_neighbor_elements_initialized) {
        size_t number_elements = indices.size();
        Faces faces;
        for (size_t i = 0; i < number_elements; i++) {
            for (int p = 0; p < dim + 1; p++) {
                int my_index = indices[i](p);
                StdVector<int>& list = vertex_to_neighbor_elements[my_index];
                auto found = std::find(std::begin(list), std::end(list), i);
                if (found == list.end())
                    list.emplace_back(i);
            }
        }
        vertex_to_neighbor_elements_initialized = true;
    }
}

/**
       Input: index is an element index.
       Input: target is a list of element indicies in the same mesh.
       Output: result is the closest one in the target to the index.
               where the distance is judged by graph connectivity.
    */
template <int dim>
void SimplexMesh<dim>::closestElement(const int index, int& result, const HashTable<int, bool>& target)
{
    if (!neighboring_initialized)
        this->initializeNeighboringElements();
    const int num_nodes = indices.size();
    bool visited[num_nodes];
    for (int i = 0; i < num_nodes; i++) {
        visited[i] = false;
    }

    std::queue<int> index_q;
    index_q.push(index);

    while (!index_q.empty()) {
        int front = index_q.front();
        visited[front] = true;
        //if front is in target, return
        if (target.get(front) != nullptr) {
            result = front;
            return;
        }

        //else, add all its neighboring index to the queue
        for (typename SpMat::InnerIterator it(connectivity, front); it; ++it) {
            ZIRAN_ASSERT(it.value() == 1);
            if (!visited[it.col()]) {
                index_q.push(it.col());
            }
        }

        index_q.pop();
    }

    result = -1;
}

template <int dim>
void SimplexMesh<dim>::findIndexBounds(int& min, int& max) const
{
    min = std::numeric_limits<int>::max();
    max = std::numeric_limits<int>::lowest();
    for (const IV& element : indices) {
        int local_min = element.minCoeff();
        int local_max = element.maxCoeff();
        if (local_min < min)
            min = local_min;
        if (local_max > max)
            max = local_max;
    }
}

template class SimplexMesh<1>;
template class SimplexMesh<2>;
template class SimplexMesh<3>;
} // namespace ZIRAN
