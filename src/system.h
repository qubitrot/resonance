#pragma once
#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <vector>
#include <unordered_map>
#include "typedefs.h"


enum ParticleType {
    PT_None,
    PT_Boson,
    PT_Fermion
};

struct Particle {
    Particle(ParticleType t) : type(t) {}

    std::string name;
    real mass;
    ParticleType type;
    int identicality; //0 = unique
};

struct InteractionV {
    enum Type {
        None,
        Gaussian,
        Harmonic
    } type;
    bool use = false;
    real v0  = 1;
    real r0  = 1;
};

class System
{
    public:
    System();
    ~System();

    //must be called after adding all particles
    //and potentials
    void init();

    void addParticle(Particle*);
    const std::vector<Particle*>& getParticles();

    void setTrappingPotential(InteractionV);
    const InteractionV& getTrappingPotential();

    void setInteractionPotential(std::string,std::string,InteractionV);
    const InteractionV& getInteraction(std::string,std::string);

    const MatrixXr& jacobiM();
    const MatrixXr& jacobiM_inv();
    const MatrixXr& lambdaM();

    const VectorXr omega(uint,uint);

    private:
    MatrixXr jacobiTransformMatrix;
    MatrixXr jacobiTM_inv;
    MatrixXr lambdaMatrix;

    std::vector<Particle*> particles;

    std::unordered_map<
        std::pair<std::string,std::string>,
        InteractionV
    > interactionPotentials;

    InteractionV trappingPotential;
};

#endif
