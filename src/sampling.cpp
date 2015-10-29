#include <ctime>
#include "sampling.h"

SamplingDistribution::SamplingDistribution(int seed)
{
    if (seed != -1) rand.seed(seed);
    else            rand.seed(std::time(nullptr));
}

SamplingDistribution::~SamplingDistribution()
{}

SD_Uniform::SD_Uniform(real mn, real mx, int seed)
    : SamplingDistribution(seed)
    , min(mn)
    , max(mx)
{}

SD_Uniform::~SD_Uniform()
{}

real SD_Uniform::operator()()
{
    real out;
    do {
        real r = (real)rand()/(real)rand.max();
        out = (max-min)*r + min;
    } while (out == 0);

    return out;
}

SD_Gaussian::SD_Gaussian(real avg, real std, real mn, real mx, int seed)
    : SamplingDistribution(seed)
    , mean(avg)
    , stdev(std)
    , min(mn)
    , max(mx)
{}

SD_Gaussian::~SD_Gaussian()
{}

real SD_Gaussian::operator()()
{
    std::normal_distribution<real> gd(mean,stdev);

    real out;
    do {
        out = gd(rand);
    } while (out <= 0 || (min != max && (out < min || out > max)));

    return out;
}

MatrixStrain::MatrixStrain(System* sys)
    : system(sys)
    , particles(sys->getParticles())
{}

MatrixStrain::~MatrixStrain()
{
    for (auto a : distributions) {
        if (a.second != nullptr) {
            delete a.second;
            a.second = nullptr;
        }
    }
}

void MatrixStrain::setDistribution(std::string& p1, std::string& p2, SamplingDistribution* sd)
{
    auto key = std::make_pair(p1,p2);
    distributions[key] = sd;
}

void MatrixStrain::computeCG(CGaussian* cg, System* sys)
{
#ifdef DEBUG_BUILD
    assert(cg->widths.rows() == cg->widths.cols());
#endif
    uint N = cg->widths.rows();
    cg->A.resize(N-1,N-1);

    for (uint k=0; k<N-1; ++k) {
        for (uint l=0; l<N-1; ++l) {
            real A_kl = 0;
            for (uint i=0; i<N; ++i) {
                for (uint j=0; j<i; ++j) {
                    VectorXr w = sys->omega(i,j);
                    A_kl += 1/(cg->widths(i,j)*cg->widths(i,j)) * w(k) * w(l);
                }
            }
            cg->A(k,l) = A_kl;
        }
    }

    cg->norm = std::pow( std::pow(twopi,N-1)/(2*cg->A).determinant() , -3./4.);
}

CGaussian MatrixStrain::genMatrix()
{
    CGaussian out;
    out.widths = genWidths();
    computeCG(&out,system);

    return out;
}

MatrixXr MatrixStrain::genWidths()
{
    uint N = particles.size();

    MatrixXr widths;
    widths.resize(N,N);

    for(uint k=0; k<N; ++k) {
        for (uint l=0; l<=k; ++l) {
            std::string p1 = particles[k]->name;
            std::string p2 = particles[k]->name;
            auto key = std::make_pair(p1,p2);

            if(distributions.find(key) == distributions.end()) {
                DEBUG_CERR("Width distribution no specfied for " << p1 << " and " << p2 << "\n");
                throw;
            }

            real w = (*distributions[key])();
            widths(k,l) = w;
            widths(l,k) = w;
        }
    }

    return widths;
}

SampleSpace::SampleSpace()
    : totalFreq(0)
{
    rand.seed(0);
}

SampleSpace::~SampleSpace()
{
    for (auto a : strains) {
        if (a.first != nullptr) {
            delete a.first;
            a.first = nullptr;
        }
    }
}

CGaussian SampleSpace::genMatrix(int s)
{
    if (s == -1) {
        uint n = chooseStrain();
        return strains[n].first->genMatrix();
    } else {
        assert( (uint)s < strains.size());
        return strains[s].first->genMatrix();
    }
}

uint SampleSpace::chooseStrain()
{
    uint r = rand() % totalFreq;
    uint n = 0;

    while (r > strains[n].second) r -= strains[n++].second;

    return n;
}

void SampleSpace::addStrain(MatrixStrain* ms, uint freq)
{
    totalFreq += freq;
    strains.push_back( std::make_pair(ms,freq) );
}
