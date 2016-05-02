#include <iostream>
#include <iomanip>
#include <ctime>
#include "sampling.h"

std::minstd_rand SamplingDistribution::rand(std::time(nullptr));

SamplingDistribution::SamplingDistribution(int seed)
{
    if (seed != -1) rand.seed(seed);
}

SamplingDistribution::~SamplingDistribution()
{}

SD_Uniform::SD_Uniform(real min, real max, int seed)
    : SamplingDistribution(seed)
    , minimum(min)
    , maximum(max)
{}

SD_Uniform::~SD_Uniform()
{}

real SD_Uniform::operator()()
{
    real out;
    do {
        real r = (real)rand()/(real)rand.max();
        out    = (maximum-minimum)*r + minimum;
    } while (out == 0);

    return out;
}

SD_Gaussian::SD_Gaussian(std::string nam, real avg, real std, real mstdf, bool has_min, real min,
                         bool has_max, real max, bool learn, uint hist_size, int seed)
    : SamplingDistribution(seed)
    , gaussian(avg,std)
    , name(nam)
    , has_minimum(has_min)
    , has_maximum(has_max)
    , minimum(min)
    , maximum(max)
    , mean(avg)
    , standard_deviation(std)
    , min_std_fac(mstdf)
    , learning(learn)
    , history_size(hist_size)
    , history_index(0)
{
    history_size = hist_size;
    history.resize(hist_size);
    for (uint i=0; i<hist_size; ++i) {
        history[i] = (*this)();
    }
}

SD_Gaussian::~SD_Gaussian()
{}

real SD_Gaussian::operator()()
{
    real out;
    do {
        out = gaussian(rand);
    } while (
           out <= 0
        || (has_minimum && out < minimum)
        || (has_maximum && out > maximum));

    return out;
}

void SD_Gaussian::learn(real value, real impact)
{
    if (!learning) return;
    (void)impact;

    history[history_index] = value;

    history_index++;
    history_index = history_index % history_size;

    mean = std::accumulate(history.begin(),history.end(),0.0)
         / history_size;

    real acc = 0;
    for (auto a : history)
        acc += a*a;

    standard_deviation = std::sqrt(acc/history_size - mean*mean);

    if (standard_deviation < mean * min_std_fac)
        standard_deviation = mean * min_std_fac;

    gaussian = std::normal_distribution<real>(mean,standard_deviation);
}

void SD_Gaussian::print_info()
{
    std::cout << std::setw(8) << name << " | "
              << "mean: " << std::setw(7) << mean
              << ", std: " << std::setw(7) << standard_deviation;
}

CG_Strain::CG_Strain()
{}

CG_Strain::~CG_Strain()
{}

CorrelatedGaussian CG_Strain::gen_widths(const std::vector<Particle>& particles)
{
    uint N = particles.size();

    Matrix<real> widths;
    widths.resize(N,N);

    for (uint k=0; k<N; ++k) {
        for (uint l=0; l<k; ++l) {
            short id1 = particles[k].id;
            short id2 = particles[l].id;
            auto  key = std::make_pair(id1,id2);

            if (distributions.find(key) == distributions.end()) {
                std::cout << "Distribution not specified between "
                          << particles[k].name << " and "
                          << particles[l].name << "\n";
                std::cout << distributions.size() << "\n";
                throw;
            }

            real w = (*distributions[key])();
            widths(k,l) = w;
            widths(l,k) = w;
        }
        widths(k,k) = 0;
    }

    CorrelatedGaussian out;

    for (uint k=0; k<N; ++k)
        out.particle_ids.push_back(particles[k].id);

    out.widths = widths;
    out.norm   = 1;
    out.strain = -1;

    return out;
}

void CG_Strain::learn(CorrelatedGaussian& cg, real impact)
{
    Matrix<real>& widths = cg.widths;

    for (uint m=0; m<widths.rows(); ++m) {
        for (uint n=0; n<m; ++n) {
            auto key = std::make_pair( cg.particle_ids[m], cg.particle_ids[n] );
            distributions[key]->learn( widths(m,n), impact );
        }
    }
}

void CG_Strain::set_distribution(short id1, short id2,
                                 const std::shared_ptr<SamplingDistribution> sd)
{
    distributions[ std::make_pair(id1,id2) ] = sd;
}

SampleSpace::SampleSpace()
    : total_frequency(0)
{
    rand.seed(0);
}

SampleSpace::~SampleSpace()
{}

uint SampleSpace::choose_strain()
{
    uint r = rand() % total_frequency;
    uint n = 0;

    while (r > strains[n].second) r -= strains[n++].second;

    return n;
}

CorrelatedGaussian SampleSpace::gen_widths(const std::vector<Particle>& particles,
                                           int strain)
{
    if (strain == -1)
        strain = choose_strain();

    CorrelatedGaussian cg = strains[strain].first.gen_widths(particles);
    cg.strain = strain;

    return cg;
}

void SampleSpace::learn(CorrelatedGaussian& cg, real impact)
{
    uint s = cg.strain;
    strains[s].first.learn(cg,impact);
}

void SampleSpace::add_strain(CG_Strain strain, uint freq)
{
    total_frequency += freq;
    strains.push_back( std::make_pair(strain,freq) );
}

void SampleSpace::print_strain_info(const std::vector<Particle>& particles, int strain)
{
    if (strain == -1) return;
    assert( strain < strains.size() && strain >= 0 );

    CG_Strain& s = strains[strain].first;

    for (auto sd : s.distributions) {
        std::string n1 = "";
        std::string n2 = "";
        for (auto p : particles) {
            if (p.id == sd.first.first)  n1 = p.name;
            if (p.id == sd.first.second) n2 = p.name;
        }
        std::cout << "        [ " << n1
                  << " - "        << n2 << " ] ";
        sd.second->print_info();
        std::cout << "\n";
    }
}
