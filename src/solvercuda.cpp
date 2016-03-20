#ifdef MAGMA

#include <set>
#include "solvercuda.h"

template<typename T>
SolverCUDA<T>::SolverCUDA(System* context)
    : SolverCPU<T>(context)
{
    magma_init();
}

template<typename T>
SolverCUDA<T>::~SolverCUDA()
{
    magma_finalize();
}

template<>
void SolverCUDA<real>::solve(Solution<real>& solution,
                             bool eigenvectors, real theta)
{
    (void)theta;

    bool need_eigenvectors = eigenvectors;

    Matrix<real> H_eigen = solution.K + solution.V;

    magmaDouble* H = convert_matrix_real(H_eigen);
    magmaDouble* O = convert_matrix_real(solution.O);

    magma_int_t n = H_eigen.rows();

    double* w;
    magma_dmalloc_cpu(&w,n);

    magma_int_t lwork;

    magma_int_t nb = magma_get_dsytrd_nb(n);
    if (need_eigenvectors) {
        lwork = std::max(2*n + n*nb, 1 + 6*n + 2*n*n);
    } else {
        lwork = 2*n + n*nb;
    }

    double* work;
    magma_dmalloc_cpu(&work,lwork);

    magma_int_t  liwork = 2*n;
    magma_int_t* iwork;
    magma_imalloc_cpu(&iwork,liwork);

    magma_int_t info = 0;

    magma_dsygvd(1,MagmaNoVec,MagmaUpper,n,H,n,O,n,w,work,lwork,iwork,liwork,&info);

    solution.eigenvalues.resize(n);
    for (int i=0; i<n; ++i) {
        solution.eigenvalues[i] = w[i];
    }

    solution.has_eigenvalues  = true;
    if (need_eigenvectors)
        solution.has_eigenvectors = true;

    magma_free_cpu(H);
    magma_free_cpu(O);
    magma_free_cpu(w);
    magma_free_cpu(work);
    magma_free_cpu(iwork);
}

template<>
void SolverCUDA<complex>::solve(Solution<complex>& solution,
                                bool eigenvectors, real theta)
{
    (void)solution;
    (void)theta;
}

template<typename T>
magmaDouble* SolverCUDA<T>::convert_matrix_real(Matrix<real>& M)
{
    magmaDouble* A = nullptr;

    magma_int_t cols = M.cols();
    magma_int_t rows = M.rows();

    magma_dmalloc_cpu(&A, rows*cols);

    for (int i=0; i<cols; ++i) {
        for (int j=0; j<rows; ++j) {
            A[i + j*rows] = M(j,i);
        }
    }

    return A;
}

template<typename T>
magmaDoubleComplex* SolverCUDA<T>::convert_matrix_complex(Matrix<complex>& M)
{
    (void)M;
}

template class SolverCUDA<real>;
template class SolverCUDA<complex>;

#endif
