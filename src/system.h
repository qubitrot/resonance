#pragma once

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

    ParticleType type;
    real mass;
    std::string name;
    short id;
    short identicality; //negative = unique
};

struct Interaction {
    enum Type {
        None,
        Gaussian,
        Harmonic,
    } type;

    bool on = false;
    real v0 = 1;
    real r0 = 1;
};

struct SymmetrizedCG
{
    std::vector<CorrelatedGaussian> funcs;
    std::vector<int> signs;
};

class System
{
public:
    System();
    ~System();

    //Must be called after setting up
    //particles and potentials.
    void init();

    void add_particle(Particle);
    void set_interaction(short id1, short id2, Interaction);
    void set_interaction(std::string n1, std::string n2, Interaction);

    const std::vector<Particle>& get_particles();
    const Interaction&           get_interaction(short id1, short id2);

    const Matrix<real>& get_jacobi_transformation();
    const Matrix<real>& get_jacobi_transformation_inverse();
    const Matrix<real>& get_lambda_matrix();

    Vector<real> omega(uint i, uint j);

    //Take a CG reference with widths defined
    //and generates everything else.
    void transform_cg(CorrelatedGaussian&);

    SymmetrizedCG symmetrize(CorrelatedGaussian&);

private:
    std::vector<Particle> particles;
    std::unordered_map<
        std::pair<short,short>,
        Interaction
    > interactions;

    Matrix<real> jacobi_transformation;
    Matrix<real> jacobi_transformation_inverse;
    Matrix<real> lambda_matrix;
};
