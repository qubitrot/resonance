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

struct Trap {
    enum Type {
        None,
        Harmonic
    } type;

    real w = 1;
};

struct Interaction {
    enum Type {
        None,
        Gaussian,
        MultiGaussian,
        PowerLaw,
        MultiPower
    } type;

    real v0   = 1;
    real pow  = 1;
    real r0sq = 1;
    real w    = 1;

    std::vector<real> mult_v0;
    std::vector<real> mult_r0sq;
    std::vector<real> mult_pow;
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
    void set_trapping_potential(Trap);

    const std::vector<Particle>& get_particles();
    const Interaction&           get_interaction(short id1, short id2);
    const Trap&                  get_trapping_potential();

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

    Trap trapping_potential;

    Matrix<real> jacobi_transformation;
    Matrix<real> jacobi_transformation_inverse;
    Matrix<real> lambda_matrix;
};
