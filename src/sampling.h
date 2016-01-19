#pragma once
#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <random>
#include <list>
#include "typedefs.h"
#include "system.h"

class SampleSpace;

class SamplingDistribution
{
public:
    SamplingDistribution(int seed, bool learn, uint hsize); //0 = seed by time
    virtual ~SamplingDistribution();

    virtual real operator()()=0;
    virtual void learn(real,real) {};

protected:
    std::minstd_rand rand;
    bool learning;

    uint hist_size;
    uint hist_index;
    std::vector<real> history;
};

class SD_Uniform : public SamplingDistribution
{
public:
    SD_Uniform(real min, real max, int seed=-1);
    ~SD_Uniform();

    virtual real operator()();

private:
    real min;
    real max;
};

class SD_Gaussian : public SamplingDistribution
{
public:
    //If min=max, then then there is no min/max
    SD_Gaussian(real avg, real std, real mn, real mx, real mstdf,
                int seed, bool learn, uint hsize);
    ~SD_Gaussian();

    virtual real operator()();
    virtual void learn(real,real);

private:
    real mean;
    real stdev;
    real min;
    real max;
    real min_std_fac;
};

class MatrixStrain
{
friend class SampleSpace;

public:
    MatrixStrain(System*);
    ~MatrixStrain();

    //Computes A and norm from the widths matrix
    static void computeCG(CGaussian*,System*);

    //Generates a new CGaussian according to distributions.
    CGaussian genMatrix();

    void setDistribution(std::string& p1, std::string& p2,
                         const std::shared_ptr<SamplingDistribution>& sd);

    void learn(CGaussian&, real impact);

private:
    std::unordered_map<
        std::pair<std::string,std::string>,
        std::shared_ptr<SamplingDistribution>
    > distributions;

    System* system;
    const std::vector<Particle*>& particles;

    MatrixXr genWidths();
};

class SampleSpace
{
public:
    SampleSpace();
    ~SampleSpace();

    //Generates a CGaussian, s chooses strain, s=-1 means semi-random
    CGaussian genMatrix(int s = -1);
    uint chooseStrain();

    void learn(CGaussian&, real impact);

    void addStrain(MatrixStrain*,uint);

private:
    std::vector<
        std::pair<MatrixStrain*,real>
    > strains;

    uint totalFreq;

    std::list< std::pair<uint,real> > learnFreqList;

    std::minstd_rand rand;
};

#endif
