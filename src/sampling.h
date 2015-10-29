#pragma once
#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <random>
#include "typedefs.h"
#include "system.h"

class SampleSpace;

class SamplingDistribution
{
public:
    SamplingDistribution(int seed=-1); //0 = seed by time
    virtual ~SamplingDistribution();

    virtual real operator()()=0;

protected:
    std::minstd_rand rand;
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
    SD_Gaussian(real avg, real std, real mn=0, real mx=0, int seed=-1);
    ~SD_Gaussian();

    virtual real operator()();

private:
    real mean;
    real stdev;
    real min;
    real max;
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

    void setDistribution(std::string& p1, std::string& p2, SamplingDistribution* sd);

private:
    std::unordered_map<
        std::pair<std::string,std::string>,
        SamplingDistribution*
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

    void addStrain(MatrixStrain*,uint);

    private:
    std::vector<
        std::pair<MatrixStrain*,uint>
    > strains;

    uint totalFreq;

    std::minstd_rand rand;
};

#endif
