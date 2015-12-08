#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include <tuple>
#include "typedefs.h"
#include "system.h"

struct SolverResults
{
    MatrixXc T;
    MatrixXc V;
    MatrixXr O;
    std::vector<complex>  eigenvalues;
    std::vector<VectorXc> eigenvectors;
};

class Solver
{
public:
    Solver(System*);
    virtual ~Solver();

    //Compute and solve for ALL the Hamiltonian elements.
    virtual SolverResults solve(const Basis&)=0;

    //Update/add row to already computed Hamiltonian/eigenvalues
    virtual SolverResults solveRow(const Basis&, SolverResults& cache, uint row)=0;

    //Compute and solve for a complex rotation. Requires unrotated cache;
    virtual SolverResults solveRot(const Basis&, real theta, SolverResults& unrot)=0;

    //Update/add rot to already computed rotated Hamiltonian/eigenvalues;
    virtual SolverResults solveRotRow(const Basis&, real theta, SolverResults& cache, uint row)=0;

    virtual real overlap(const CGaussian&, const CGaussian&)=0;

protected:
    System*  system;
};

class CpuSolver : public Solver
{
public:
    CpuSolver(System*);
    ~CpuSolver();

    SolverResults solve(const Basis&);
    SolverResults solveRow(const Basis&, SolverResults& cache, uint row);
    SolverResults solveRot(const Basis&, real theta, SolverResults& unrot);
    SolverResults solveRotRow(const Basis&, real theta, SolverResults& cache, uint row);

    real overlap(const CGaussian&, const CGaussian&);

private:
    complex kinetic  (const CGaussian&, const CGaussian&, real);
    complex gaussianV(real v0, real r0, real theta, real over, real c_ij);

    real genc_ij(const CGaussian& A, const CGaussian& B, uint i, uint j);

    SolverResults computeHermition(MatrixXc& T, MatrixXc& V, MatrixXr& O);
    SolverResults computeQZ(MatrixXc& T, MatrixXc& V, MatrixXr& O);

    std::tuple<Basis,std::vector<int>,uint> symmetrize(const CGaussian&);
};

#endif
