#pragma once

#ifdef MAGMA

#include "solver.h"
#include "magma.h"

typedef double magmaDouble;

template<typename T>
class SolverCUDA : public SolverCPU<T>
{
public:
    SolverCUDA(System* context);
    ~SolverCUDA();

    virtual void solve(Solution<T>&, bool eigenvectors=0, real theta=0);

protected:
    magmaDouble*        convert_matrix_real(Matrix<real>&);
    magmaDoubleComplex* convert_matrix_complex(Matrix<complex>&);
};

#endif
