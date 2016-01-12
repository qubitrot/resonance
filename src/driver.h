#pragma once
#ifndef DRIVER_H
#define DRIVER_H

#include <map>
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
    Basis generateBasis(uint size, bool rot=false, real start=0, real end=0, uint steps=0);

    void sweepAngle(real start, real end, uint steps);
    void updateSweep(uint row);

    void readBasis(std::string file, uint n=0, bool append=false);
    void writeBasis(std::string file);
    void writeConvergenceData(std::string file);
    void writeSweepData(std::string file);
    void printEnergies(uint n);

    int  targetState;
    real targetEnergy;
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
    std::map<real,SolverResults> sweepData;
    std::tuple<real,real,uint> sweepMetaData; //start,end,steps

    static void findBestAddition(std::pair<CGaussian,complex>* out,Driver*,Basis trails,
                                 SolverResults* bcache,uint target,real singularityLimit);
};

#endif
