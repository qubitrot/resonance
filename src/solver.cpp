#include <set>
#include <algorithm>
#include "solver.h"
#include "sampling.h"

Solver::Solver(System* sys)
    : system(sys)
{}

Solver::~Solver()
{}

CpuSolver::CpuSolver(System* sys)
    : Solver(sys)
{}

CpuSolver::~CpuSolver()
{}

SolverResults CpuSolver::solve(const Basis& basis, real theta)
{
    const std::vector<Particle*>& particles = system->getParticles();

    uint size = basis.size();
    uint N    = particles.size();

    MatrixXc H;     H.resize(size,size);
    MatrixXr O;     O.resize(size,size);

    for (uint m=0; m<size; ++m) {
        for (uint n=0; n<=m; ++n) {
            H(m,n) = 0;
            O(m,n) = 0;

            auto              symFunc = symmetrize(basis[m]);
            Basis&            A_sym   = std::get<0>(symFunc);
            std::vector<int>& signs   = std::get<1>(symFunc);
            uint              nperm   = std::get<2>(symFunc);

            for (uint k=0; k<A_sym.size(); ++k) {
                real ol = overlap(A_sym[k],basis[n]);

                O(m,n) += signs[k]*nperm * ol;
                H(m,n) += complex(signs[k]*nperm) * kinetic(A_sym[k],basis[n],theta,ol);

                for (uint i=0; i<N; ++ i) {
                    for (uint j=0; j<i; ++j) {
                        std::string p1 = particles[i]->name;
                        std::string p2 = particles[j]->name;

                        real c_ij = genc_ij(A_sym[k],basis[n],i,j);

                        const InteractionV& V = system->getInteraction(p1,p2);
                        switch (V.type) {
                            case InteractionV::Gaussian:
                                H(m,n) += complex(signs[k]*nperm)
                                        * gaussianV(V.v0,V.r0,theta,ol,c_ij);
                                break;

                            case InteractionV::Harmonic:
                                break;
                            case InteractionV::None:
                                break;
                        }
                    }
                }
            }
            O(n,m) = O(m,n);
            H(n,m) = H(m,n);
        }
    }

    if (theta == 0) return computeHermition(H,O);
    else            return compute(H,O);
}

SolverResults CpuSolver::solveRow(const Basis& basis, real theta, SolverResults& cache, uint row)
{
    const std::vector<Particle*>& particles = system->getParticles();

    uint N    = particles.size();
    uint size = basis.size();

#ifdef DEBUG_BUILD
    assert(cache.H.rows() == size ||
           cache.H.rows() == size -1);
    assert(row < size);

    if (cache.H.rows() == size-1) {
        assert(row+1 == size);
    }
#endif

    MatrixXc H = cache.H;
    MatrixXr O = cache.O;
    H.conservativeResize(size,size);
    O.conservativeResize(size,size);

    uint n = row;
    for (uint m=0; m<size; ++m) {
        H(m,n) = 0;
        O(m,n) = 0;

        auto              symFunc = symmetrize(basis[m]);
        Basis&            A_sym   = std::get<0>(symFunc);
        std::vector<int>& signs   = std::get<1>(symFunc);
        uint              nperm   = std::get<2>(symFunc);

        for (uint k=0; k<A_sym.size(); ++k) {
            real ol = overlap(A_sym[k],basis[n]);

            O(m,n) += signs[k]*nperm * ol;
            H(m,n) += complex(signs[k]*nperm) * kinetic(A_sym[k],basis[n],theta,ol);

            for (uint i=0; i<N; ++ i) {
                for (uint j=0; j<i; ++j) {
                    std::string p1 = particles[i]->name;
                    std::string p2 = particles[j]->name;

                    real c_ij = genc_ij(A_sym[k],basis[n],i,j);

                    const InteractionV& V = system->getInteraction(p1,p2);
                    switch (V.type) {
                        case InteractionV::Gaussian:
                            H(m,n) += complex(signs[k]*nperm)
                                    * gaussianV(V.v0,V.r0,theta,ol,c_ij);
                            break;

                        case InteractionV::Harmonic:
                            break;
                        case InteractionV::None:
                            break;
                    }
                }
            }
        }

        O(n,m) = O(m,n);
        H(n,m) = H(m,n);
    }

    if (theta == 0) return computeHermition(H,O);
    else            return compute(H,O);
}

SolverResults CpuSolver::compute(MatrixXc& H, MatrixXr& O)
{
    Eigen::ComplexEigenSolver<MatrixXc> eigenSolver;
    eigenSolver.compute(H*O.inverse().cast<complex>());

    VectorXc eigenvals = eigenSolver.eigenvalues().cast<complex>();

    SolverResults out;
    out.H = H;
    out.O = O;

    out.eigenvalues.resize(eigenvals.rows());
    for (int k=0; k<eigenvals.rows(); ++k) {
        out.eigenvalues[k] = eigenvals(k);
    }

    struct {
        bool operator()(complex a, complex b) {
            return a.real() < b.real();
        }
    } customLess;
    std::sort(out.eigenvalues.begin(),out.eigenvalues.end(),customLess);

    return out;
}

SolverResults CpuSolver::computeHermition(MatrixXc& H, MatrixXr& O)
{
    Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXc> eigenSolver;
    eigenSolver.compute(H,O.cast<complex>());

    VectorXc eigenvals = eigenSolver.eigenvalues().cast<complex>();

    SolverResults out;
    out.H = H;
    out.O = O;

    out.eigenvalues.resize(eigenvals.rows());
    for (int k=0; k<eigenvals.rows(); ++k) {
        out.eigenvalues[k] = eigenvals(k);
    }

    struct {
        bool operator()(complex a, complex b) {
            return a.real() < b.real();
        }
    } customLess;
    std::sort(out.eigenvalues.begin(),out.eigenvalues.end(),customLess);

    return out;
}

/*SolverResults CpuSolver::MAGMAcomputeHermition(MatrixXc& H, MatrixXc& O)
{
    double* w;
    double* A;
    double* B;
    double* work;

    magma_int_t n    = H.rows();
    magma_int_t code = magma_dsygvd(1,MagmaNoVec,MagmaUpper,n,A,n,B,n,work,);
}*/

real CpuSolver::overlap(const CGaussian& A, const CGaussian& B)
{
    register constexpr real twopi = 2*pi;
    register uint n = A.A.rows();
    register real q = std::pow(twopi,n)/(A.A+B.A).determinant();

    return A.norm*B.norm * q * std::sqrt(q);
}

complex CpuSolver::kinetic(const CGaussian& A, const CGaussian& B, real theta, real over)
{
    register constexpr real k = 3./2. * hbar;

    MatrixXr C = A.A.selfadjointView<Eigen::Lower>()
               * (A.A+B.A).inverse()
               * B.A.selfadjointView<Eigen::Lower>() * system->lambdaM();

    complex  v = complex(k*C.trace()*over,0);

    return std::exp(complex(0,-2)*theta)*v;
}

complex CpuSolver::gaussianV(real v0, real r0, real theta, real over, real c_ij)
{
    constexpr real k = pi*sqrt(pi);

    complex a = complex(1./2.*c_ij,0) + std::exp(complex(0,2*theta)) / complex(2*r0*r0,0);
    complex integral = complex(v0*k,0) * std::pow(a,-3./2.);

    return std::pow(c_ij/(2*pi), 3./2.) * over * integral;
}

real CpuSolver::genc_ij(const CGaussian& A, const CGaussian& B, uint i, uint j)
{
    MatrixXr C = (A.A+B.A).inverse();
    VectorXr w_ij = system->omega(i,j);
    real c_ij = 1./( w_ij.transpose() * C * w_ij );

    return c_ij;
}

std::tuple<Basis,std::vector<int>,uint> CpuSolver::symmetrize(const CGaussian& A)
{
    struct Group {
        std::vector<uint> indices;
        ParticleType type;
    };

    struct Permutation {
        std::vector<uint> indices;
        int sign;
    };

    std::vector<Particle*> particles = system->getParticles();
    uint N = particles.size();

    //Identify all the groups
    std::vector<Group> groups;

    std::set<int> assigned;
    for (uint k=0; k<N; ++k) {
        //If this element hasn't been assigned to a group yet, create one and add it
        if (assigned.find(k) != assigned.end()) continue;

        Group group;
        group.indices.push_back(k);
        group.type = particles[k]->type;

        //find all the other particles with same identicality and add it to group
        Particle* p1 = particles[k];
        for (uint l=k+1; l<N; ++l) {
            Particle* p2 = particles[l];

            if (p1->identicality == p2->identicality) {
                group.indices.push_back(l);
                assigned.insert(l);
            }
        }

        groups.push_back(group);
    }

    //Permute the particles in each group
    Permutation defaultPerm;
    defaultPerm.sign = 1;
    defaultPerm.indices.resize(N);
    for (uint k=0; k<N; ++k)
        defaultPerm.indices[k]=k;

    std::vector<Permutation> permutations;
    permutations.push_back(defaultPerm);

    for (Group group : groups) {
        int sign = 1;
        Group gperm = group;
        while (std::next_permutation(gperm.indices.begin(),gperm.indices.end())) {
            if (group.type == ParticleType::PT_Fermion) sign *= -1;

            Permutation perm;
            perm.sign = sign;
            perm.indices.resize(N);

            uint gindex = 0;
            for (uint n=0; n<N; ++n) {
                if (n == group.indices[gindex]) {
                    perm.indices[n] = gperm.indices[gindex];
                    gindex++;
                } else {
                    perm.indices[n] = n;
                }
            }

            permutations.push_back(perm);
        }
    }

    //Operate permutations on matrices
    Basis symA;
    std::vector<int> signs;

    for (auto perm : permutations) {
        MatrixXr P;
        P.resize(N,N);
        for (uint i=0; i<N; ++i) {
            for (uint j=0; j<=i; ++j) {
                P(i,j) = A.widths(perm.indices[i],perm.indices[j]);
                P(j,i) = P(i,j);
            }
        }

        CGaussian cg;
        cg.widths = P;
        MatrixStrain::computeCG(&cg,system);

        symA.push_back(cg);
        signs.push_back(perm.sign);
    }

    /*std::cout << A.A << "\n";
    std::cout << "PERMUTED:\n\n";
    for (auto B : symA) {
        std::cout << B.A << "\n\n";
    }*/

    return std::make_tuple(symA,signs,permutations.size());
}

