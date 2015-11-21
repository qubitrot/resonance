#pragma once
#ifndef DRIVER_H
#define DRIVER_H

#include "typedefs.h"
#include "system.h"
#include "sampling.h"
#include "solver.h"

class Driver
{
public:
    Driver(System*,Solver*,SampleSpace*);
    ~Driver();

    Basis generateTrials(uint size);
    Basis generateBasis(uint size);

    void sweepAngle(real start, real end, uint steps);

    void readBasis(std::string file, uint n=0, bool append=false);
    void writeBasis(std::string file);
    void writeConvergenceData(std::string file);

    uint targetState;
    uint trialSize;
    uint numThreads;
    real singularityLimit;
    bool forceDiversity;

private:
    System*      system;
    Solver*      solver;
    SampleSpace* sampleSpace;

    Basis         basis;
    SolverResults basisCache;

    std::vector<complex> convergenceData;

    static void findBestAddition(std::pair<CGaussian,complex>* out,Driver*,Basis trails,
                                 SolverResults* bcache,uint target,real singularityLimit);
};

#endif
