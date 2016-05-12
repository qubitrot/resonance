#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <stack>
#include <thread>
#include <mutex>
#include "driver.h"
#include "json/json.h"
#include "profiler/profiler.h"

Driver::Driver(System* sys, SampleSpace* ss)
    : target_state(0)
    , target_energy(-10)
    , targeting_energy(false)
    , trial_size(100)
    , singularity_limit(1e-10)
    , threads(1)
    , system(sys)
    , sample_space(ss)
{}

Driver::~Driver()
{}

ConvergenceData Driver::expand_basis(Basis& basis, uint size)
{
    Solution<real> empty_cache;
    return expand_basis(basis,empty_cache,size);
}

ConvergenceData Driver::expand_basis(Basis& basis, Solution<real>& cache, uint size)
{
    PROFILE();

    SolverCPU<real> solver(system);

    ConvergenceData out;
    out.offset = basis.size();

    uint fail_svd_count = 0;

    for (uint s=0; s<size; ++s) {
        std::mutex mtx;
        std::stack<
            std::pair<CorrelatedGaussian,
                      Solution<real>>
        > trial_stack;

        uint target = target_state;
        if (basis.size() < target_state)
            target = basis.size();

        auto trial_thread = [&]() {
            Basis trials = generate_trials(trial_size);

            for (uint i=0; i<trial_size; ++i) {
                CorrelatedGaussian trial = trials[i];

                Basis trial_basis = basis;
                trial_basis.push_back(trial);

                Solution<real> trial_solution = solver.solve(trial_basis,cache);

                mtx.lock();
                if (trial_stack.empty() ||
                    trial_solution.eigenvalues[target] <
                    trial_stack.top().second.eigenvalues[target]) {
                        trial_stack.push( std::make_pair(trial,trial_solution) );
                }
                //std::cout << " " << trial_solution.eigenvalues[target] << " ";
                //std::flush(std::cout);
                mtx.unlock();
            }
        };

        if (threads == 1) {
            trial_thread();
        } else {
            std::vector<std::thread> thread_vec;
            for (uint i=0; i<threads; ++i) {
                thread_vec.push_back(std::thread(trial_thread));
            }
            for (uint i=0; i<threads; ++i) {
                thread_vec[i].join();
            }
        }

        uint num_candidates = trial_stack.size();
        int  chosen_strain  = -1;
        real lowest_ev;
        uint svd_rejections = 0;

        //Make sure there is not too much linear dependance
        bool everything_fails_SVD = true;
        while (!trial_stack.empty()) {
            Solution<real> trial_solution = trial_stack.top().second;

            lowest_ev = lowest_eigenval(trial_solution.O);
            if (lowest_ev > singularity_limit) {
                basis.push_back(trial_stack.top().first);
                cache = trial_solution;
                sample_space->learn(trial_stack.top().first,0);
                chosen_strain = trial_stack.top().first.strain;
                everything_fails_SVD = false;
                break;
            }
            trial_stack.pop();
            svd_rejections++;
        }

        std::cout << "------------------------------------------------------------------------\n";
        std::cout << "Addition of basis function #" << basis.size() << "\n\n";
        std::cout << "      Target E:    ";
        if (targeting_energy) std::cout << std::setw(10) << target_energy;
        else                  std::cout << std::setw(10) << "NA";
        std::cout << "      Trials: " << trial_size*threads
                  << ", Candidates: " << num_candidates << "\n";
        std::cout << "      Target St.   " << std::setw(10) << target
                  << "      SVD rejections: " << svd_rejections << "\n";

        for (uint i=0; i<14; ++i) std::cout << "\n";

        if (everything_fails_SVD) {
            for (uint i=0; i<2; ++i) std::cout << "\x1b[A";
            std::cout << "!     All candidates failed SVD.";
            fail_svd_count++;
            s--;

            if (fail_svd_count > 1) {
                std::cout << "   Stepping back...\n\n";
                if (!basis.empty())           basis.pop_back();
                if (!out.eigenvalues.empty()) out.eigenvalues.pop_back();
                cache = solver.solve(basis);
                s--;
            } else {
                std::cout << "\n\n";
            }

        } else {
            solver.solve(cache,true); //get eigenvectors;
            out.eigenvalues.push_back( cache.eigenvalues );

            if (targeting_energy && target == target_state &&
                    cache.eigenvalues[target] < target_energy)
                target_state++;

            fail_svd_count = 0;

            for (uint i=0; i<14; ++i) std::cout << "\x1b[A";

            std::cout << "      Eigenenergy: " << std::setw(10) << std::setprecision(5)
                                               << cache.eigenvalues[target]
                      << "      Singularity: " << lowest_ev << "\n\n";
            real csp = 0;
            for (auto b : basis) {
                if (b.strain == chosen_strain)
                    csp++;
            }
            csp *= 100.0/basis.size();

            std::cout << "      Chose strain " << chosen_strain << ", info:        "
                      << "Prevalence: "  << std::setw(4) << std::setprecision(3) << csp << "%\n";
            sample_space->print_strain_info(system->get_particles(), chosen_strain);

            Eigen::IOFormat CleanFmt(4,0, ", ", "\n", "        [", "]");

            std::cout << "\n      Widths of new function:\n"
                      << basis[basis.size()-1].widths.format(CleanFmt)
                      << "\n\n";
        }

    }

    return out;
}

real Driver::lowest_eigenval(Matrix<real>& O)
{
    PROFILE();

    /*Eigen::JacobiSVD<Matrix<real>> svd(O);

    real lowest_sv = svd.singularValues()(0);
    for (uint i=0; i<svd.singularValues().rows(); ++i) {
        if (svd.singularValues()(i) < lowest_sv) {
            lowest_sv = svd.singularValues()(i);
        }
    }*/

    Eigen::SelfAdjointEigenSolver<Matrix<real>> eigen_solver(O);

    real lowest = eigen_solver.eigenvalues()(0);
    for (uint i=0; i<eigen_solver.eigenvalues().rows(); ++i) {
        if (eigen_solver.eigenvalues()(i) < lowest)
            lowest = eigen_solver.eigenvalues()(i);
    }

    return lowest;
}

Basis Driver::generate_trials(uint n)
{
    PROFILE();

    Basis out;
    const std::vector<Particle>& particles = system->get_particles();

    for (uint i=0; i<n; ++i) {
        CorrelatedGaussian cg = sample_space->gen_widths(particles);
        system->transform_cg(cg);

        out.push_back(cg);
    }

    return out;
}

SweepData Driver::sweep_basis(Basis& basis, real start, real end, uint steps)
{
    Solution<complex> empty_cache;
    return sweep_basis(basis,empty_cache,start,end,steps);
}

SweepData Driver::sweep_basis(Basis& basis, Solution<complex>& cache,
                              real start, real end, uint steps)
{
    PROFILE();

    SolverCPU<complex> solver(system);

    std::mutex mtx;
    SweepData out;

    real step_size = (end-start)/steps;

    auto sweep_thread = [&](real theta) {
        Solution<complex> sol = solver.solve(basis,cache,false,theta);

        sol.eigenvectors.clear();
        sol.K.resize(0,0);
        sol.V.resize(0,0);
        sol.O.resize(0,0);

        mtx.lock();
        out.sweep_vec.push_back( std::make_pair(theta,sol) );
        mtx.unlock();

    };

    if (threads == 1) {
        for (uint i=0; i<steps; ++i) {
            real theta = start + i*step_size;
            std::cout << "theta = " << theta << " | ";
            std::flush(std::cout);
            sweep_thread(theta);
            std::cout << "done.\n";
        }
    } else {
        for (uint i=0; i<steps;) {
            std::vector<std::thread> thread_vec;
            std::cout << "theta = ";
            for (uint j=0; j<threads && i<steps;) {
                real theta = start + i*step_size;
                thread_vec.push_back(
                        std::thread(sweep_thread,theta)
                );
                std::cout << theta << ", ";
                j++;
                i++;
            }
            std::cout << "\n";
            for (uint j=0; j<thread_vec.size(); ++j) {
                thread_vec[j].join();
            }
        }
    }

    out.start     = start;
    out.end       = end;
    out.steps     = steps;
    out.step_size = step_size;

    return out;
}

PairDistribution Driver::pair_distribution(Basis& basis, Solution<real>& sol, uint i, uint j,
                                           real start, real num_boxes)
{
    real r = 1;

    for (uint m=0; m<basis.size(); ++m) {
    }
}

Basis Driver::read_basis(std::string file, uint n)
{
    PROFILE();

    Json::Value  root;
    Json::Reader reader;

    std::ifstream basis_file(file, std::ifstream::binary);
    if (!basis_file.good()) {
        throw;
    }

    basis_file >> root;

    Basis basis;

    Json::Value set = root["basis"];
    if (n == 0) n = set.size();

    for (uint k=0; k<set.size() && k<n; ++k) {
        Json::Value entry = set[k];

        uint n = std::sqrt(entry["A"].size());
        Matrix<real> A;
        A.resize(n,n);

        for (uint i=0; i<n; ++i) {
            for (uint j=0; j<n; ++j) {
                A(i,j) = entry["A"][n*i + j].asDouble();
            }
        }

        uint m = std::sqrt(entry["widths"].size());
        Matrix<real> widths;
        widths.resize(m,m);

        for (uint i=0; i<m; ++i) {
            for (uint j=0; j<m; ++j) {
                widths(i,j) = entry["widths"][m*i + j].asDouble();
            }
        }

        assert(m == n+1);

        CorrelatedGaussian cg;
        cg.trans  = A;
        cg.widths = widths;
        cg.norm   = entry["norm"].asDouble();

        basis.push_back(cg);
    }

    return basis;
}

void Driver::write_basis(Basis& basis, std::string file)
{
    PROFILE();

    Json::Value root;

    Json::Value metadata(Json::objectValue);
    metadata["size"] = Json::Value::UInt(basis.size());
    root["metadata"] = metadata;

    Json::Value set(Json::arrayValue);
    for (auto cg : basis) {
        Json::Value entry(Json::objectValue);

        Matrix<real> trans = cg.trans;
        Json::Value A(Json::arrayValue);
        for (int i=0; i<trans.rows(); ++i) {
            for (int j=0; j<trans.cols(); ++j) {
                A.append(Json::Value( trans(i,j) ));
            }
        }

        Matrix<real> widths = cg.widths;
        Json::Value w_matrix(Json::arrayValue);
        for (int i=0; i<widths.rows(); ++i) {
            for (int j=0; j<widths.cols(); ++j) {
                w_matrix.append(Json::Value( widths(i,j) ));
            }
        }

        entry["A"]      = A;
        entry["widths"] = w_matrix;
        entry["norm"]   = cg.norm;
        entry["strain"] = cg.strain;

        set.append(entry);
    }

    root["basis"] = set;

    std::ofstream basis_file;
    basis_file.open(file, std::ofstream::out);
    basis_file << root;
    basis_file.close();
}

void Driver::write_convergence(ConvergenceData& cd, std::string file,bool append)
{
    PROFILE();

    std::ofstream datafile;

    if (append) datafile.open(file, std::ofstream::app);
    else        datafile.open(file, std::ofstream::out);

    for (uint i=0; i<cd.eigenvalues.size(); ++i) {
        datafile << cd.offset + i;
        for (auto ev : cd.eigenvalues[i]) {
            datafile << "\t" << ev;
        }
        datafile << "\n";
    }

    datafile.close();
}

void Driver::write_sweep(SweepData& sd, std::string file, bool append)
{
    PROFILE();

    std::ofstream datafile;

    if (append) datafile.open(file, std::ofstream::app);
    else        datafile.open(file, std::ofstream::out);

    for (uint i=0; i<sd.steps; ++i) {
        std::pair<real,Solution<complex>> sol = sd.sweep_vec[i];
        datafile << sol.first << "\t";
        for (auto ev : sol.second.eigenvalues) {
            datafile << ev.real() << " " << ev.imag() << "\t";
        }
        datafile << "\n";
    }

    datafile.close();
}
