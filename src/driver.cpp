#include <thread>
#include <fstream>
#include <iomanip>
#include "driver.h"
#include "json/json.h"

Driver::Driver(System* sys, Solver* sol, SampleSpace* ss)
    : targetState(0)
    , trialSize(0)
    , numThreads(1)
    , overlapLimit(0.98)
    , forceDiversity(false)
    , system(sys)
    , solver(sol)
    , sampleSpace(ss)
{}

Driver::~Driver()
{}

Basis Driver::generateTrials(uint n)
{
    Basis trials;

    for (uint s=0; s<n; ++s) {
        CGaussian cg;

        int strainN = -1;
        if (forceDiversity) strainN = sampleSpace->chooseStrain();

        bool overlapProblem;
        do {
            cg = sampleSpace->genMatrix(strainN);
            overlapProblem = false;
            for (unsigned k=0; k<basis.size(); ++k) {
                if (solver->overlap(cg,basis[k]) > overlapLimit) {
                    overlapProblem = true;
                    break;
                }
            }
        } while (overlapProblem);

        trials.push_back(cg);
    }

    return trials;
}

Basis Driver::generateBasis(uint size)
{
    if (basisCache.eigenvalues.size() != basis.size()) {
        basisCache = solver->solve(basis,0);
    }

    for (uint s=0; s<size; ++s) {
        std::cout << basis.size() << " | ";

        std::vector<std::pair<CGaussian,complex>*> candidates;
        std::vector<SolverResults*> caches;
        std::vector<std::thread> threads;

        std::cout << "target: " << targetState << "    ";

        for (uint i=0; i<numThreads; ++i) {
            candidates.push_back(new std::pair<CGaussian,complex>);
            caches.push_back(new SolverResults);
            *caches[i] = basisCache;
            threads.push_back(
                    std::thread(Driver::findBestAddition,candidates[i],this,
                                generateTrials(trialSize),caches[i],targetState,0)
            );
        }

        auto it = threads.begin();
        for (; it != threads.end(); ++it) {
            it->join();
        }

        CGaussian best = candidates[0]->first;
        complex   bev  = candidates[0]->second;
        basisCache = *caches[0];
        for (uint i=0; i<candidates.size(); ++i) {
            auto c = candidates[i];
            if (c->second.real() < bev.real()) {
                best = c->first;
                bev  = c->second;
                basisCache = *caches[i];
            }
            delete c;
            delete caches[i];
        }

        basis.push_back(best);
        convergenceData.push_back(bev);

        if (basis.size() > targetState) {
            std::cout << "E = " << std::setprecision(18) << bev << "\n";
        } else {
            std::cout << "\n";
        }

        //std::cout << "____\n";
        //for (auto b : basis) std::cout << "A:\n" << b.A << "\n\n";
    }

    return basis;
}

void Driver::findBestAddition(std::pair<CGaussian,complex>* out, Driver* driver, Basis trials,
                              SolverResults* bcache, uint target, real theta)
{
    CGaussian best;
    complex lowestEV = complex(0,0);

    SolverResults cache = *bcache;

    if (target >= driver->basis.size())
        target  = driver->basis.size();

    for (auto cg : trials) {
        Basis test = driver->basis;
        test.push_back(cg);

        cache = driver->solver->solveRow(test,theta,cache,test.size()-1);

        std::vector<complex> ev = cache.eigenvalues;
        if (lowestEV == complex(0,0) || ev[target].real() < lowestEV.real()) {
            out->first  = cg;
            out->second = ev[target];
            *bcache     = cache;
        }
    }
}

void Driver::sweepAngle(uint steps, real stepsize)
{
    std::vector<SolverResults> vsr;
    vsr.resize(steps);

    for (uint i=0; i<steps; ++i) {
        std::vector<std::thread> threads;

        for(uint n=0; n<numThreads; ++n) {
            real theta = stepsize*i;

            if (i >= steps) break;

            std::cout << i << " Solve angle " << theta << "\n";
            threads.push_back(threadify_member(&Solver::solve,solver,&vsr[i],basis,theta));

            i++;
        }
        i--;

        auto it = threads.begin();
        for (; it != threads.end(); ++it) {
            it->join();
        }
    }

    std::ofstream sweepFile("sweep.dat");

    for (uint i=0; i<steps; ++i) {
        sweepFile << stepsize*i;
        for (auto a : vsr[i].eigenvalues)
            sweepFile << "\t" << a.real() << " " << a.imag();
        sweepFile << "\n";
    }
}

void Driver::writeBasis(std::string file)
{
    Json::Value root;

    Json::Value metadata(Json::objectValue);
    metadata["size"] = Json::Value::UInt(basis.size());
    root["metadata"] = metadata;

    Json::Value set(Json::arrayValue);
    for (auto cg : basis) {
        Json::Value entry(Json::objectValue);

        MatrixXr A = cg.A;
        Json::Value aMatrix(Json::arrayValue);
        for (int i=0; i<A.rows(); ++i) {
            for (int j=0; j<A.cols(); ++j) {
                aMatrix.append(Json::Value( A(i,j) ));
            }
        }

        MatrixXr widths = cg.widths;
        Json::Value wMatrix(Json::arrayValue);
        for (int i=0; i<widths.rows(); ++i) {
            for (int j=0; j<widths.cols(); ++j) {
                wMatrix.append(Json::Value( widths(i,j) ));
            }
        }

        entry["A"]      = aMatrix;
        entry["widths"] = wMatrix;
        entry["norm"]   = cg.norm;

        set.append(entry);
    }

    root["basis"] = set;

    std::ofstream basisFile;
    basisFile.open(file, std::ofstream::out);
    basisFile << root;
    basisFile.flush();
}

void Driver::readBasis(std::string file)
{
    Json::Value  root;
    Json::Reader reader;

    std::ifstream basisFile(file, std::ifstream::binary);
    if (!basisFile.good()) {
        return;
    }

    basisFile >> root;

    Basis newBasis;

    Json::Value set = root["basis"];
    for (uint k=0; k<set.size(); ++k) {
        Json::Value entry = set[k];

        uint n = std::sqrt(entry["A"].size());
        MatrixXr A;
        A.resize(n,n);

        for (uint i=0; i<n; ++i) {
            for (uint j=0; j<n; ++j) {
                A(i,j) = entry["A"][n*i + j].asDouble();
            }
        }

        uint m = std::sqrt(entry["widths"].size());
        MatrixXr widths;
        widths.resize(m,m);

        for (uint i=0; i<m; ++i) {
            for (uint j=0; j<m; ++j) {
                widths(i,j) = entry["widths"][m*i + j].asDouble();
            }
        }

        assert(m == n+1);

        CGaussian cg;
        cg.A      = A;
        cg.widths = widths;
        cg.norm   = entry["norm"].asDouble();

        newBasis.push_back(cg);
    }

    basis = newBasis;
}

void Driver::writeConvergenceData(std::string file)
{
    std::ofstream datafile;
    datafile.open(file, std::ofstream::out);

    for (uint i=0; i<convergenceData.size(); ++i) {
        datafile << i << "\t" << convergenceData[i].real()
                 << " "       << convergenceData[i].imag() <<  "\n";
    }

    datafile.flush();
}
