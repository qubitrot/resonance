#include <set>
#include "solver.h"
#include "profiler/profiler.h"

typedef struct {
    double real;
    double imag;
} fortranComplex16;

extern "C" void zggev_(const char* JOBVL, const char* JOBVR, const int* N,
                       fortranComplex16* A, const int* LDA, fortranComplex16* B,
                       const int* LDB, fortranComplex16* ALPHA, fortranComplex16* BETA,
                       fortranComplex16* VL, const int* LDVL, fortranComplex16* VR,
                       const int* LDVR, fortranComplex16* WORK, const int* LWORK,
                       double* RWORK, int* INFO);

template<typename T>
SolverCPU<T>::SolverCPU(System* context)
    : Solver<T>(context)
{}

template<typename T>
SolverCPU<T>::~SolverCPU()
{}

template<typename T>
Solution<T> SolverCPU<T>::compute(Basis& basis, real theta)
{
    Solution <T> empty_cache;
    return compute(basis, empty_cache, theta);
}

template<typename T>
Solution<T> SolverCPU<T>::compute(Basis& basis, Solution<T>& cache, real theta)
{
    PROFILE();

    auto particles = this->system_context->get_particles();

    uint size = basis.size();
    uint N    = particles.size();

    //TODO: do NOT assume that cache is correct;
    Matrix<T>    K(cache.K);  K.conservativeResize(size,size);
    Matrix<T>    V(cache.V);  V.conservativeResize(size,size);
    Matrix<real> O(cache.O);  O.conservativeResize(size,size);

    uint initial_index = cache.K.rows();
    if (theta != cache.theta)
        initial_index = 0;

    for (uint m=initial_index; m<size; ++m) {

        SymmetrizedCG symcg = this->system_context->symmetrize(basis[m]);
        auto funcs   = symcg.funcs;
        uint nperm   = funcs.size();

        for (uint n=0; n<=m; ++n) {
            K(m,n) = 0;
            V(m,n) = 0;
            O(m,n) = 0;

            for (uint k=0; k<funcs.size(); ++k) {
                real ol = overlap(funcs[k],basis[n]);

                real prefac = symcg.signs[k] * nperm;

                O(m,n) += prefac * ol;
                K(m,n) += prefac * kinetic(funcs[k],basis[n],ol);

                //particle interactions
                for (uint i=0; i<N; ++i) {
                    for(uint j=0; j<i; ++j) {
                        short id1 = particles[i].id;
                        short id2 = particles[j].id;

                        const Interaction& inter =
                            this->system_context->get_interaction(id1,id2);

                        real cij;

                        switch (inter.type) {
                            case Interaction::Gaussian:
                                cij = c_ij(funcs[k],basis[n],i,j);
                                V(m,n) += prefac * gaussian_v(inter.v0,inter.r0sq,ol,cij,theta);
                                break;
                            case Interaction::MultiGaussian:
                                cij = c_ij(funcs[k],basis[n],i,j);
                                for (uint a=0; a<inter.mult_v0.size(); ++a) {
                                     real v0   = inter.mult_v0[a];
                                     real r0sq = inter.mult_r0sq[a];
                                     V(m,n) += prefac * gaussian_v(v0,r0sq,ol,cij,theta);
                                }
                                break;
                            default:
                                break;
                        }
                    }
                }
            }

            rotate_if_T_complex(K(m,n),-1*theta);

            K(n,m) = K(m,n);
            V(n,m) = V(m,n);
            O(n,m) = O(m,n);

        }
    }

    Solution<T> solution;
    solution.K = K;
    solution.V = V;
    solution.O = O;
    solution.theta = theta;
    solution.eigenvalues  = cache.eigenvalues;
    solution.eigenvectors = cache.eigenvectors;
    solution.has_eigenvalues  = false;
    solution.has_eigenvectors = false;

    return solution;
}

template<typename T>
Solution<T> SolverCPU<T>::solve(Basis& basis, bool eigenvectors, real theta)
{
    Solution<T> empty_cache;
    Solution<T> sol = compute(basis,empty_cache,theta);
    solve(sol,eigenvectors);
    return sol;
}

template<typename T>
Solution<T> SolverCPU<T>::solve(Basis& basis, Solution<T>& cache,
                                bool eigenvectors, real theta)
{
    Solution<T> sol = compute(basis,cache,theta);
    solve(sol,eigenvectors);
    return sol;
}

template<>
void SolverCPU<real>::solve_bisection(Solution<real>& solution,
                                      uint max_iterations, real tolorance)
{
    PROFILE();

    uint size = solution.K.rows();

    real norm = solution.O(size-1,size-1);
    for (uint i=0; i<size-1; ++i) {
        real ol = 0;
        for (uint j=0; j<size-1; ++j) {
            ol += solution.eigenvectors[i](j)
                * solution.O(j,size-1);
        }
        norm -= std::abs(ol*ol);
    }

    Matrix<real> H = Matrix<real>::Zero(size,size);

    for (uint i=0; i<solution.eigenvalues.size(); ++i) {
        H(i,i) = solution.eigenvalues[i];
    }

    for (uint i=0; i<size-1; ++i) {
        for (uint j=0; j<size-1; ++j) {
            H(i,size-1) += solution.eigenvectors[i](j)
                        *( (solution.K(j,size-1) + solution.V(j,size-1))
                         - (solution.eigenvalues[i] * solution.O(j,size-1)));
        }
        H(i,size-1) *= 1./std::sqrt(norm);
        H(size-1,i) = H(i,size-1);
    }

    H(size-1,size-1) = solution.K(size-1,size-1) + solution.V(size-1,size-1);
    for (uint i=0; i<size-1; ++i) {
        real fac1 = 0;
        real fac2 = 0;
        for (uint j=0; j<size-1; ++j) {
            fac1 += solution.eigenvectors[i](j) * solution.O(size-1,j);
            fac2 += solution.eigenvectors[i](j) * (
                        solution.K(j,size-1) + solution.V(j,size-1)
                    );
        }
        H(size-1,size-1) -= 2*fac1*fac2;
        H(size-1,size-1) += solution.eigenvalues[i] * std::abs(fac1*fac1);
    }
    H(size-1,size-1) *= 1./norm;

    //Root finding time.
    auto polynomial = [&H,&size,&solution](real E) {
        real out = H(size-1,size-1) - E;
        for (uint i=0; i<size-1; ++i) {
            out -= std::abs( H(i,size-1)*H(i,size-1) ) / (solution.eigenvalues[i] - E);
        }
        return out;
    };

    std::vector<real> new_eigenvalues;

    for (uint k=0; k<solution.eigenvalues.size(); ++k) {
        real E_a;
        real E_b  = solution.eigenvalues[k];
             E_b += 0.000001;

        if (k>0) {
            E_a  = solution.eigenvalues[k-1];
            E_a += 0.000001;
        } else {
            E_a = solution.eigenvalues[k]-100;
        }

        real midpoint = (E_a + E_b)/2;

        for (uint n=0; n<max_iterations; ++n) {
            midpoint = (E_a + E_b)/2;
            if ((E_b - E_a)/2.0 < tolorance) {
                break;
            }

            real p_a = polynomial(E_a);
            real p_m = polynomial(midpoint);

            if ( ((p_a > 0) && (p_m > 0)) || ((p_a < 0) && (p_m < 0)) ) {
                E_a = midpoint;
            } else {
                E_b = midpoint;
            }
        }
        new_eigenvalues.push_back(midpoint);
    }

    solution.eigenvalues = new_eigenvalues;
    solution.has_eigenvalues  = true;
}

template<>
void SolverCPU<complex>::solve_bisection(Solution<complex>& solution,
                                         uint max_iterations, real tolorance)
{
    //This should never be called;
    assert(false);
    (void)solution;
    (void)max_iterations;
    (void)tolorance;
}

template<>
void SolverCPU<real>::solve_full(Solution<real>& solution, bool eigenvectors)
{
    PROFILE();

    (void)eigenvectors;

    Matrix<real> H = solution.K + solution.V;

    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix<real>> eigen_solver;
    eigen_solver.compute(H,solution.O);

    Vector<real> eigenvals = eigen_solver.eigenvalues();
    Matrix<real> eigenvecs = eigen_solver.eigenvectors();

    //TODO: Do sorting better
    std::vector<std::pair<real,int>> sort_vec;
    for (uint k=0; k<eigenvals.rows(); ++k) {
        sort_vec.push_back( std::make_pair(eigenvals(k),k) );
    }

    struct {
        bool operator()(std::pair<real,int> a, std::pair<real,int> b) {
            return a.first < b.first;
        }
    } custom_less;
    std::sort(sort_vec.begin(),sort_vec.end(),custom_less);

    solution.eigenvalues.resize( eigenvals.rows() );
    solution.eigenvectors.resize( eigenvecs.rows() );
    for (int k=0; k<eigenvals.rows(); ++k) {
        solution.eigenvalues[k]  = sort_vec[k].first;
        solution.eigenvectors[k] = eigenvecs.col( sort_vec[k].second );
    }

    solution.has_eigenvalues  = true;
    solution.has_eigenvectors = true;
}

template<>
void SolverCPU<complex>::solve_full(Solution<complex>& solution, bool eigenvectors)
{
    PROFILE();

    (void)eigenvectors;

    Matrix<complex> H = solution.K + solution.V;

    int n     = H.rows();
    int ldvl  = n;
    int ldvr  = n;
    int lwork = 2*n;
    int info;

    fortranComplex16* A     = new fortranComplex16[n*n];
    fortranComplex16* B     = new fortranComplex16[n*n];
    fortranComplex16* alpha = new fortranComplex16[n];
    fortranComplex16* beta  = new fortranComplex16[n];
    fortranComplex16* work  = new fortranComplex16[std::max(1,lwork)];
    double*           rwork = new double[8*n];

    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            A[n*i+j].real = H.data()[n*i+j].real();
            A[n*i+j].imag = H.data()[n*i+j].imag();
            B[n*i+j].real = solution.O.data()[n*i+j];
            B[n*i+j].imag = 0;
        }
    }

    zggev_("N","N",&n,A,&n,B,&n,alpha,beta,nullptr,&ldvl,nullptr,&ldvr,
            work,&lwork,rwork,&info);

    std::vector<complex> eigenvalues;
    for (int i=0; i<n; ++i) {
        eigenvalues.push_back(complex(alpha[i].real,alpha[i].imag)
                             /complex(beta[i].real,beta[i].imag));
    }

    struct {
        bool operator()(complex a, complex b) {
            return a.real() < b.real();
        }
    } customLess;
    std::sort(eigenvalues.begin(),eigenvalues.end(),customLess);

    delete[] rwork;
    delete[] work;
    delete[] beta;
    delete[] alpha;
    delete[] B;
    delete[] A;

    solution.eigenvalues = eigenvalues;
    solution.has_eigenvalues = true;
}

template<>
void SolverCPU<real>::solve(Solution<real>& solution, bool eigenvectors)
{
    PROFILE();

    if (!eigenvectors && solution.eigenvectors.size() > 0 &&
        (uint)solution.K.rows() == solution.eigenvectors.size() +1 &&
        (uint)solution.K.rows() == solution.eigenvalues.size()  +1 )
    {
        return solve_bisection(solution,10e6,10e-6);
    } else {
        return solve_full(solution,eigenvectors);
    }
}

template<>
void SolverCPU<complex>::solve(Solution<complex>& solution, bool eigenvectors)
{
    PROFILE();

    return solve_full(solution,eigenvectors);
}

template<typename T>
inline
real SolverCPU<T>::overlap(CorrelatedGaussian& A, CorrelatedGaussian& B)
{
    PROFILE();

    constexpr real twopi = 2*pi;
    uint n = A.trans.rows();
    real q = std::pow(twopi,n)/(A.trans + B.trans).determinant();

    return A.norm * B.norm * q * std::sqrt(q);
}

template<typename T>
inline
T SolverCPU<T>::kinetic(CorrelatedGaussian& A, CorrelatedGaussian& B, real over)
{
    PROFILE();

    constexpr real k = 3./2. * hbar;

    Matrix<real> C = A.trans * (A.trans+B.trans).inverse()
                   * B.trans * this->system_context->get_lambda_matrix();

    return T( k * C.trace() * over );
}

template<typename T>
inline
T SolverCPU<T>::gaussian_v(real v0, real r0sq, real over, real cij, real theta)
{
    PROFILE();

    (void)theta;
    constexpr real k = 5.568327996831707845284817982118; //pi * sqrt(pi);

    T r = 1;
    rotate_if_T_complex(r,theta);

    T a        = T( 1./2. * cij + r/(2*r0sq) );
    T integral = T( v0*k ) * std::pow(a,-3./2.);

    real x = cij/(2*pi);

    return x * std::sqrt(x) * over * integral;
}

template<typename T>
inline
real SolverCPU<T>::c_ij(CorrelatedGaussian& A, CorrelatedGaussian& B, uint i, uint j)
{
    PROFILE();

    Matrix<real> S = A.trans + B.trans;
    Matrix<real> C = S.llt().solve(Matrix<real>::Identity(S.rows(),S.cols()));
    Vector<real> w = this->system_context->omega(i,j);

    real a = (w.transpose() * C * w);

    return 1./a;
}

template<>
inline
void SolverCPU<real>::rotate_if_T_complex(real& value, real theta)
{}

template<>
inline
void SolverCPU<complex>::rotate_if_T_complex(complex& value, real theta)
{
    value *= std::exp(complex(0,2*theta));
}

template class SolverCPU<real>;
template class SolverCPU<complex>;
