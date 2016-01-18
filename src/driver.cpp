#include <thread>
#include <fstream>
#include <iomanip>
#include "driver.h"
#include "json/json.h"

Driver::Driver(System* sys, Solver* sol, SampleSpace* ss)
    : targetState(0)
    , targetEnergy(1111)
    , trialSize(0)
    , numThreads(1)
    , singularityLimit(5e-14)
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
    int strainN = -1;
    if (forceDiversity) strainN = sampleSpace->chooseStrain();

    for (uint s=0; s<n; ++s) {
        CGaussian cg;
        cg = sampleSpace->genMatrix(strainN);
        trials.push_back(cg);
    }

    return trials;
}

Basis Driver::generateBasis(uint size, bool rot, real start, real end, uint steps)
{
    if (basisCache.eigenvalues.size() != basis.size()) {
        basisCache = solver->solve(basis);
    }

    for (uint s=0; s<size; ++s) {
        std::cout << basis.size() << " | ";

        std::vector<std::pair<CGaussian,complex>*> candidates;
        std::vector<SolverResults*> caches;
        std::vector<std::thread> threads;

        while (targetEnergy != 1111 &&
            basisCache.eigenvalues.size() > targetState &&
            basisCache.eigenvalues[targetState].real() < targetEnergy) {
                targetState++;
        }

        std::cout << "target: " << targetState << "  ";
        std::flush(std::cout);

        for (uint i=0; i<numThreads; ++i) {
            candidates.push_back(new std::pair<CGaussian,complex>);
            caches.push_back(new SolverResults);
            *caches[i] = basisCache;
            threads.push_back(
                    std::thread(Driver::findBestAddition,candidates[i],this,generateTrials(trialSize),
                                caches[i],targetState,singularityLimit)
            );
        }

        auto it = threads.begin();
        for (; it != threads.end(); ++it) {
            it->join();
        }

        bool redo = true; //to accomidate the possibility that all trials
                          //exhibit too much linear dependance.

        CGaussian best = candidates[0]->first;
        complex   bev  = candidates[0]->second;
        basisCache = *caches[0];
        for (uint i=0; i<candidates.size(); ++i) {
            auto c = candidates[i];
            if (c->first.A.rows() != 0 && c->second.real() <= bev.real()) {
                redo = false;
                best = c->first;
                bev  = c->second;
                basisCache = *caches[i];
            }
            delete c;
            delete caches[i];
        }

        if (!redo) {
            basis.push_back(best);
            convergenceData.push_back(bev);

            std::cout << "strain: " << best.strain << "  ";

            if (basis.size() > targetState + 1) {
                std::cout << "E = " << std::setprecision(18) << bev << "\n";
            } else {
                std::cout << "\n";
            }

            if (rot) {
                if (sweepMetaData == std::make_tuple(start,end,steps)) {
                    updateSweep(basis.size()-1);
                } else {
                    sweepAngle(start,end,steps);
                }
            }

            sampleSpace->learn(best,0);

        } else {
            std::cout << "REDO: Too much linear dependance.\n";
            s--;
        }

    }

    return basis;
}

void Driver::findBestAddition(std::pair<CGaussian,complex>* out, Driver* driver, Basis trials,
                              SolverResults* bcache, uint target, real singularityLimit)
{
    CGaussian best;
    complex lowestEV = complex(0,0);

    SolverResults cache = *bcache;

    if (target >= driver->basis.size())
        target  = driver->basis.size();

    out->first.A.resize(0,0); //"uninitialized" used for redo check

    std::stack< std::tuple<CGaussian,complex,SolverResults> > stack;

    for (auto cg : trials) {
        Basis test = driver->basis;
        test.push_back(cg);

        cache = driver->solver->solveRow(test,cache,test.size()-1);

        std::vector<complex> ev = cache.eigenvalues;
        if (lowestEV == complex(0,0) || ev[target].real() < lowestEV.real()) {
            stack.push(
                    std::make_tuple(cg,ev[target],cache)
            );
            lowestEV = ev[target];
        }
    }

    while (!stack.empty()) {
        auto candidate = stack.top();
        MatrixXr& O    = std::get<2>(candidate).O;

        //Make sure there's not too much linear dependance
        Eigen::JacobiSVD<MatrixXr> svd(O);
        real lowestSV = svd.singularValues()(0);
        for (uint i=0; i<svd.singularValues().rows(); ++i) {
            if (svd.singularValues()(i) < lowestSV) lowestSV = svd.singularValues()(i);
        }

        if (lowestSV > singularityLimit) {
            out->first  = std::get<0>(candidate);
            out->second = std::get<1>(candidate);
            *bcache     = std::get<2>(candidate);

            break;
        }

        stack.pop();
    }
}

void Driver::sweepAngle(real start, real end, uint steps)
{
    assert (basis.size() > 0);

    //sweepMetaData = std::make_tuple(start,end,steps);
    //sweepData.clear();

    real stepsize = (end-start)/(real)steps;

    if (basisCache.eigenvalues.size() != basis.size()) {
        basisCache = solver->solve(basis);
    }

    if (numThreads == 1) {
        for (uint i=0; i<steps; ++i) {
            real theta = start + stepsize*i;
            std::cout << i << " Solve angle " << theta << "\n";

            sweepData[theta] = solver->solveRot(basis,theta,basisCache);
        }
    } else {
        for (uint i=0; i<steps; ++i) {
            std::vector<std::thread> threads;

            for(uint n=0; n<numThreads; ++n) {
                real theta = start + stepsize*i;

                if (i >= steps) break;

                std::cout << i << " Solve angle " << theta << "\n";
                threads.push_back(threadify(
                                  &Solver::solveRot,&sweepData[theta],solver,basis,theta,basisCache));
                i++;
            }
            i--;

            auto it = threads.begin();
            for (; it != threads.end(); ++it) {
                it->join();
            }
        }
    }
}

void Driver::updateSweep(uint row)
{
    assert (basis.size() > 0);

    auto i = sweepData.begin();

    if (numThreads == 1) {
        for (; i != sweepData.end(); ++i) {
            real theta = i->first;
            SolverResults& cache = i->second;

            sweepData[theta] = solver->solveRotRow(basis,theta,cache,row);
        }
    } else {
        for (; i != sweepData.end(); ++i) {
            std::vector<std::thread> threads;

            for(uint n=0; n<numThreads; ++n) {
                if (i == sweepData.end()) break;

                real theta = i->first;
                SolverResults& cache = i->second;

                threads.push_back(threadify(
                                  &Solver::solveRotRow,&sweepData[theta],solver,basis,theta,cache,row));
                ++i;
            }
            --i;

            auto it = threads.begin();
            for (; it != threads.end(); ++it) {
                it->join();
            }
        }
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
    basisFile.close();
}

void Driver::readBasis(std::string file, uint n, bool append)
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
    if (n == 0) n = set.size();

    for (uint k=0; k<set.size() && k<n; ++k) {
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

    if (append) basis.insert(basis.end(), newBasis.begin(), newBasis.end());
    else        basis = newBasis;
}

void Driver::writeConvergenceData(std::string file)
{
    std::ofstream datafile;
    datafile.open(file, std::ofstream::out);

    for (uint i=0; i<convergenceData.size(); ++i) {
        datafile << i+basis.size()-convergenceData.size() << "\t" << convergenceData[i].real()
                 << " " << convergenceData[i].imag() <<  "\n";
    }

    datafile.close();
}

void Driver::writeSweepData(std::string file)
{
    std::ofstream datafile;
    datafile.open(file, std::ofstream::out);

    for (auto& kv : sweepData) {
        datafile << kv.first;
        for (auto a : kv.second.eigenvalues)
            datafile << "\t" << a.real() << " " << a.imag();
        datafile << "\n";
    }

    datafile.close();
}

void Driver::printEnergies(uint n)
{
    if (basisCache.eigenvalues.size() != basis.size()) {
        basisCache = solver->solve(basis);
    }

    for (uint i=0; i<basisCache.eigenvalues.size() && i<n; ++i) {
        std::cout << "E" << i << " = ";
        std::cout << std::setprecision(18) << basisCache.eigenvalues[i];
        std::cout << "\n";
    }
}
