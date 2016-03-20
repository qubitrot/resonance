#pragma once

#include <string>
#include "typedefs.h"
#include "system.h"
#include "solver.h"
#include "solvercuda.h"
#include "sampling.h"

struct ConvergenceData {
    uint offset;
    std::vector<
        std::vector<real>
    > eigenvalues;
};

struct SweepData {
    real start;
    real end;
    real steps;
    real step_size;
    std::vector<
        std::pair<real,Solution<complex>> //theta, solution
    > sweep_vec;
};

class Driver
{
public:
    Driver(System*,SampleSpace*);
    ~Driver();

    ConvergenceData expand_basis(Basis& basis, uint size);
    ConvergenceData expand_basis(Basis& basis, Solution<real>& cache, uint size);

    SweepData sweep_basis(Basis& basis, real start, real end, uint steps);
    SweepData sweep_basis(Basis& basis, Solution<complex>& cache,
                          real start, real end, uint steps);
    SweepData sweep_basis(Basis& basis, real start, real end, uint steps,
                          std::vector<uint> sizes);
    SweepData sweep_basis(Basis& basis, Solution<complex>& cache, real start,
                          real end, uint steps, std::vector<uint> sizes);

    Basis generate_trials(uint n);

    static Basis read_basis(std::string file, uint n=0);
    static void  write_basis(Basis& basis, std::string file);
    static void  write_convergence(ConvergenceData& cd, std::string file,
                                   bool append = false);
    static void  write_sweep(SweepData& sd, std::string file,
                             bool append = false);

    uint target_state;
    real target_energy;
    bool targeting_energy;
    uint trial_size;
    real singularity_limit;
    uint threads;

private:
    System*      system;
    SampleSpace* sample_space;

    bool check_SVD(Matrix<real>& O);
};
