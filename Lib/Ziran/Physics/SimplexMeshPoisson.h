#ifndef SIMPLEX_MESH_POISSON_H
#define SIMPLEX_MESH_POISSON_H

#include <Eigen/Sparse>
#include <Ziran/Math/Geometry/SimplexMesh.h>
#include <Ziran/Math/Linear/KrylovSolvers.h>
#include <Ziran/Math/Linear/Minres.h>
#include <algorithm>

/**
This is a class that solves poisson equation -\nabla u = f.
Discretization: mat* x =b.
**/

namespace ZIRAN {

template <class T, int dim>
class SimplexMeshPoisson {
public:
    typedef Vector<T, dim> TV;
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, Eigen::Dynamic> VEC;

    SimplexMesh<dim>& tet_mesh;
    StdVector<TV>& X;
    StdVector<T> element_measure;
    StdVector<T> b_scale;
    StdVector<Eigen::Matrix<T, dim, dim + 1>> grad_N;
    T unit_simplex_measure;
    Eigen::SparseMatrix<T, Eigen::RowMajor> mat;
    int n; //The total number of dofs
    VEC& x;
    VEC& b;
    StdVector<int>& dirichlet;
    StdVector<T>& dirichlet_value;

    SimplexMeshPoisson(SimplexMesh<dim>& mesh,
        StdVector<TV>& X,
        VEC& x,
        VEC& b,
        StdVector<int>& dirichlet_in,
        StdVector<T>& dirichlet_value_in)
        : tet_mesh(mesh)
        , X(X)
        , element_measure(mesh.numberElements())
        , b_scale(mesh.numberElements())
        , grad_N(mesh.numberElements())
        , n(X.size())
        , x(x)
        , b(b)
        , dirichlet(dirichlet_in)
        , dirichlet_value(dirichlet_value_in)
    {
        ZIRAN_ASSERT(n > 0);
        T d_factorial = (T)1;
        for (int i = 0; i < dim; i++)
            d_factorial *= (T)(i + 1);
        T unit_simplex_measure = (T)1 / d_factorial;

        //set up Dm_inverse and interpolating function derivatives (wrt rest configuration)
        TM dm;
        T scale = (T)1 / (T)((dim + 1) * (dim + 2));
        for (size_t e = 0; e < tet_mesh.numberElements(); e++) {
            dS(e, dm);
            element_measure[e] = unit_simplex_measure * dm.determinant();
            b_scale[e] = element_measure[e] * scale;
            TM dm_inv_trans = dm.inverse().transpose();
            grad_N[e] << -dm_inv_trans.rowwise().sum(), dm_inv_trans;
        }

        //build matrix
        mat.resize(n, n);
        mat.reserve(Eigen::VectorXi::Constant(n, 27).transpose());
        for (size_t e = 0; e < tet_mesh.numberElements(); e++) {
            for (int i = 0; i <= dim; i++) {
                for (int j = 0; j <= dim; j++) {
                    mat.coeffRef(tet_mesh.index(e, i), tet_mesh.index(e, j)) += grad_N[e].col(i).dot(grad_N[e].col(j)) * element_measure[e];
                }
            }
        }
        mat.makeCompressed();

        initialize_x();
    }

    void initialize_x()
    {
        x.setZero();
        for (size_t i = 0; i < dirichlet.size(); i++)
            x(dirichlet[i]) = dirichlet_value[i];
    }

    void dS(const size_t element, TM& result)
    {
        for (int i = 0; i < dim; i++)
            result.col(i) = X[tet_mesh.index(element, i + 1)] - X[tet_mesh.index(element, 0)];
    }

    //set b based on the analytic rhs
    void setRHS(const VEC& rhs)
    {
        b.resize(X.size());
        b.setZero();
        for (size_t e = 0; e < tet_mesh.numberElements(); e++) {
            for (int i = 0; i <= dim; i++) {
                for (int j = 0; j <= dim; j++) {
                    T this_scale = b_scale[e] * ((i == j) + 1);
                    b(tet_mesh.index(e, i)) += this_scale * rhs(tet_mesh.index(e, j));
                }
            }
        }
    }

    int solve()
    {
        EigenSparseKrylovMatrix<T, VEC> eska(mat);

        DirichletDofProjection<VEC> projection_function(dirichlet);
        eska.setProjection(projection_function);

        JacobiPreconditioner<VEC> preconditioner(mat);
        eska.setPreconditioner(preconditioner);

        Minres<T, EigenSparseKrylovMatrix<T, VEC>, VEC> minres(10000);
        return minres.solve(eska, x, b, false);
    }

    //after solving the poission for \u(x), get \grad \u(x) for each tetrahedron
    void findGradient(StdVector<Vector<T, dim>>& gradPhi)
    {
        gradPhi.resize(tet_mesh.numberElements());
        for (size_t e = 0; e < tet_mesh.numberElements(); e++) {
            gradPhi[e].setZero();
            for (int i = 0; i <= dim; i++)
                gradPhi[e] += x(tet_mesh.index(e, i)) * grad_N[e].col(i);
        }
    }

    //after solving the poission for \u(x), get normalized \grad \u(x) for each tetrahedron
    void findNormalizedGradient(StdVector<Vector<T, dim>>& gradPhi)
    {
        gradPhi.resize(tet_mesh.numberElements());
        for (size_t e = 0; e < tet_mesh.numberElements(); e++) {
            gradPhi[e].setZero();
            for (int i = 0; i <= dim; i++)
                gradPhi[e] += x(tet_mesh.index(e, i)) * grad_N[e].col(i);

            if (gradPhi[e].norm() < 16 * std::numeric_limits<T>::epsilon())
                gradPhi[e].setZero();
            else
                gradPhi[e].normalize();
        }
    }
};
} // namespace ZIRAN

#endif
