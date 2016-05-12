#pragma once

#include <unordered_map>
#include <memory>
#include <random>
#include "typedefs.h"
#include "system.h"

class SampleSpace;

class SamplingDistribution
{
public:
    SamplingDistribution(int seed=-1); //seed by time if neg
    virtual ~SamplingDistribution();

    virtual real operator()()=0;

    virtual void learn(real value, real impact) {
        (void)value;
        (void)impact;
    };

    virtual void print_info() {};

protected:
    static std::minstd_rand rand;
};

class SD_Uniform : public SamplingDistribution
{
public:
    SD_Uniform(real min, real max, int seed=-1);
    ~SD_Uniform();

    virtual real operator()();

private:
    real minimum;
    real maximum;
};

class SD_Gaussian : public SamplingDistribution
{
public:
    SD_Gaussian(std::string name, real avg, real std, real mstdf, bool has_min, real min,
                bool has_max, real max, bool learn, uint hist_size, int seed=-1);
    ~SD_Gaussian();

    virtual real operator()();
    virtual void learn(real value, real impact);
    virtual void print_info();

private:
    std::normal_distribution<real> gaussian;

    std::string name;

    bool has_minimum;
    bool has_maximum;
    real minimum;
    real maximum;
    real mean;
    real standard_deviation;
    real min_std_fac;

    bool learning;
    uint history_size;
    uint history_index;
    std::vector<real> history;
};

class CG_Strain
{
friend class SampleSpace;
public:
    CG_Strain();
    ~CG_Strain();

    CorrelatedGaussian gen_widths(const std::vector<Particle>& particles);
    void learn(CorrelatedGaussian& cg, real impact);
    void set_distribution(short id1, short id2,
                          const std::shared_ptr<SamplingDistribution> sd);

private:
    std::unordered_map<
        std::pair<short,short>,
        std::shared_ptr<SamplingDistribution>
    > distributions;
};

class SampleSpace
{
public:
    SampleSpace();
    ~SampleSpace();

    uint choose_strain();
    CorrelatedGaussian gen_widths(const std::vector<Particle>& particles,
                                  int strain = -1);

    void learn(CorrelatedGaussian& cg, real impact);

    void add_strain(CG_Strain, uint freq);

    void print_strain_info(const std::vector<Particle>& particles, int strain=-1);

private:
    std::vector<
        std::pair<CG_Strain,real>
    > strains;

    uint total_frequency;

    std::minstd_rand rand;
};




