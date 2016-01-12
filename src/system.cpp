#include <numeric>
#include "system.h"

System::System()
{}

System::~System()
{
    for (auto p : particles) {
        if (p != nullptr) delete p;
        p = nullptr; //be safe
    }
}

void System::init()
{
    int N = particles.size();

    //Generate the Jacobi transformation matrix
    std::vector<real> masses;
    for (auto p : particles) {
        masses.push_back(p->mass);
    }

    MatrixXr U;
    U.resize(N,N);

    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
                 if (j-i >  1) U(i,j) = 0;
            else if (j-i == 1) U(i,j) = -1;
            else               U(i,j) = masses[j] /
                                        std::accumulate(masses.begin(),masses.begin()+i+1,0.d);
        }
    }

    jacobiTransformMatrix = U;
    jacobiTM_inv          = U.inverse();

    //Generate lambda matrix
    MatrixXr L;
    L.resize(N-1,N-1);

    for (int i=0; i<N-1; ++i) {
        for (int j=0; j<N-1; ++j) {
            real L_ij = 0;
            for (int k=0; k<N; ++k) {
                L_ij += U(i,k)*U(j,k) / masses[k];
            }
            L(i,j) = L_ij;
        }
    }

    lambdaMatrix = L;
}

void System::addParticle(Particle* p)
{
    particles.push_back(p);

    //set default interactions (none)
    for (auto p2 : particles) {
        setInteractionPotential(p->name, p2->name, InteractionV());
    }
}

const std::vector<Particle*>& System::getParticles()
{
    return particles;
}

void System::setTrappingPotential(InteractionV V)
{
    trappingPotential = V;
}

const InteractionV& System::getTrappingPotential()
{
    return trappingPotential;
}

void System::setInteractionPotential(std::string p1, std::string p2, InteractionV V)
{
    auto key = std::make_pair(p1,p2);
    interactionPotentials[key] = V;
}

const InteractionV& System::getInteraction(std::string p1, std::string p2)
{
   auto key = std::make_pair(p1,p2);

#ifdef DEBUG_BUILD
    if (interactionPotentials.find(key) != interactionPotentials.end()) {
        return interactionPotentials[key];
    } else {
        std::cerr << "Request for nonexistent interaction potential.\n";
        throw;
    }
#endif

    return interactionPotentials[key];
}

const MatrixXr& System::jacobiM()
{
    return jacobiTransformMatrix;
}

const MatrixXr& System::jacobiM_inv()
{
    return jacobiTM_inv;
}

const MatrixXr& System::lambdaM()
{
    return lambdaMatrix;
}

const VectorXr System::omega(uint i, uint j)
{
    uint N = particles.size();

    VectorXr w;
    w.resize(N-1);

    for (uint k=0; k<N-1; ++k) {
        w(k) = jacobiTM_inv(i,k) - jacobiTM_inv(j,k);
    }

    return w;
}


