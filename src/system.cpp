#include <set>
#include <numeric>
#include "system.h"
#include "profiler/profiler.h"

System::System()
{}

System::~System()
{}

void System::init()
{
    int N = particles.size();

    //Generate the Jacobi transformation matrix;
    std::vector<real> masses;
    for (auto p : particles)
        masses.push_back(p.mass);

    Matrix<real> U;
    U.resize(N,N);

    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
                 if (j-i >  1) U(i,j) =  0;
            else if (j-i == 1) U(i,j) = -1;
            else
                 U(i,j) = masses[j] /
                          std::accumulate(masses.begin(),masses.begin()+i+1,0.0);
        }
    }

    jacobi_transformation         = U;
    jacobi_transformation_inverse = U.inverse();

    //Generate lambda matrix
    Matrix<real> L;
    L.resize(N-1,N-1);

    for (int i=0; i<N-1; ++i) {
        for (int j=0; j<N-1; ++j) {
            L(i,j) = 0;
            for (int k=0; k<N; ++k)
                L(i,j) += U(i,k)*U(j,k) / masses[k];
        }
    }

    lambda_matrix = L;
}

void System::add_particle(Particle p)
{
    particles.push_back(p);
    for (auto p2 : particles) {
        Interaction nope; nope.type = Interaction::Type::None;
        interactions[std::make_pair(p.id,p2.id)] = nope;
    }
}

void System::set_interaction(short id1, short id2, Interaction v)
{
    interactions[std::make_pair(id1,id2)] = v;
}

void System::set_interaction(std::string n1, std::string n2, Interaction v)
{
    short id1 = -999;
    short id2 = -999;

    for (auto p : particles) {
        if (p.name == n1) id1 = p.id;
        if (p.name == n2) id2 = p.id;
    }

    if (id1 == -999 || id2 == -999) throw;

    set_interaction(id1,id2,v);
}

const std::vector<Particle>& System::get_particles()
{
    return particles;
}

const Interaction& System::get_interaction(short id1, short id2)
{
    return interactions[std::make_pair(id1,id2)];
}

const Matrix<real>& System::get_jacobi_transformation()
{
    return jacobi_transformation;
}

const Matrix<real>& System::get_jacobi_transformation_inverse()
{
    return jacobi_transformation_inverse;
}

const Matrix<real>& System::get_lambda_matrix()
{
    return lambda_matrix;
}

Vector<real> System::omega(uint i, uint j)
{
    PROFILE();

    uint N = particles.size();

    Vector<real> w;
    w.resize(N-1);

    for (uint k=0; k<N-1; ++k) {
        w(k) = jacobi_transformation_inverse(i,k)
             - jacobi_transformation_inverse(j,k);
    }

    return w;
}

void System::transform_cg(CorrelatedGaussian& cg)
{
    PROFILE();

    uint N = cg.widths.rows();

    cg.trans.resize(N-1,N-1);

    for (uint k=0; k<N-1; ++k) {
        for (uint l=0; l<N-1; ++l) {
            cg.trans(k,l) = 0;
            for (uint i=0; i<N; ++i) {
                for (uint j=0; j<i; ++j) {
                    Vector<real> w = omega(i,j);
                    cg.trans(k,l) += 1./(cg.widths(i,j)*cg.widths(i,j)) * w(k) * w(l);
                }
            }
        }
    }

    cg.norm = std::pow( std::pow(2*pi,N-1)/(2*cg.trans).determinant(), -3./4.);
}

SymmetrizedCG System::symmetrize(CorrelatedGaussian& A)
{
    PROFILE();

    struct Group {
        std::vector<uint> indices;
        ParticleType type;
    };

    struct Permutation {
        std::vector<uint> indices;
        int sign;
    };

    uint N = particles.size();

    //Identify all the groups
    std::vector<Group> groups;

    std::set<int> assigned;
    for (uint k=0; k<N; ++k) {
        //If this element hasn't been assigned to a group yet, create one and add it
        if (assigned.find(k) != assigned.end()) continue;

        Group group;
        group.indices.push_back(k);
        group.type = particles[k].type;

        //find all the other particles with same identicality and add it to group
        Particle& p1 = particles[k];
        for (uint l=k+1; l<N; ++l) {
            Particle& p2 = particles[l];

            if (p1.identicality == p2.identicality) {
                group.indices.push_back(l);
                assigned.insert(l);
            }
        }

        groups.push_back(group);
    }

    //Permute the particles in each group
    Permutation default_perm;
    default_perm.sign = 1;
    default_perm.indices.resize(N);
    for (uint k=0; k<N; ++k)
        default_perm.indices[k]=k;

    std::vector<Permutation> permutations;
    permutations.push_back(default_perm);

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
    Basis sym_A;
    std::vector<int> signs;

    for (auto perm : permutations) {
        Matrix<real> P;
        P.resize(N,N);
        for (uint i=0; i<N; ++i) {
            for (uint j=0; j<=i; ++j) {
                P(i,j) = A.widths(perm.indices[i],perm.indices[j]);
                P(j,i) = P(i,j);
            }
        }

        CorrelatedGaussian cg;
        cg.widths = P;
        transform_cg(cg);

        sym_A.push_back(cg);
        signs.push_back(perm.sign);
    }

    SymmetrizedCG out;

    out.funcs = sym_A;
    out.signs = signs;

    return out;
}
