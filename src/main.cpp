#include <iostream>
#include <fstream>
#include <ctime>
#include "system.h"
#include "solver.h"
#include "sampling.h"
#include "driver.h"
#include "json/json.h"
#include "profiler/profiler.h"
#include "profiler/htmlwriter.h"

void init(std::string, System*, Driver*, SampleSpace*);

int main(int argc, char* argv[])
{
#ifndef NDEBUG
    Profiler::setLogFile("out/preformance/report.html");
    Profiler::setLogFormat(new HTMLWriter);
    PROFILE();
#endif

    std::string out_dir     = "";
    std::string config_file = "";
    std::string config_name = "";

    if (argc < 2) {
        std::cout << "Nope.\n";
        return 1;
    } else {
        config_file = argv[1];
        if (argc > 2) {
            out_dir = argv[2];
        } else {
            time_t t = time(0);
            struct tm* now = localtime(&t);

            std::string date;
            date = std::to_string(now->tm_year+1900) + std::string("-")
                 + std::to_string(now->tm_mon+1)     + std::string("-")
                 + std::to_string(now->tm_mday)      + std::string("-")
                 + std::to_string(now->tm_hour);

            int lastdot   = config_file.find_last_of(".");
            int lastslash = config_file.find_last_of("/");

            config_name = config_file.substr(lastslash+1,lastdot-lastslash-1);
            out_dir     = "out/" + date + "_" + config_name;
        }
    }

    std::string mkdir = std::string("mkdir -p ") + out_dir;
    std::system(mkdir.c_str());
    std::system(("cp " + config_file + " " + out_dir
                 + "/" + config_name + ".json").c_str());

    System* system            = new System();
    SampleSpace* sample_space = new SampleSpace();
    Driver* driver            = new Driver(system,sample_space);

    init(config_file,system,driver,sample_space);

    Json::Value  root;
    Json::Reader reader;
    std::ifstream config_doc(config_file, std::ifstream::binary);
    config_doc >> root;
    const Json::Value config_queue = root["queue"];

    SolverCPU<real> solver(system);
    Basis           basis;
    Solution<real>  cache;

    for (uint i=0; i<config_queue.size(); ++i) {
        const Json::Value item = config_queue[i];
        std::string action = item[0].asString();

        if (action == "read") {
            std::string file = item[1].asString();
            uint num = item[2].asInt();
            basis    = driver->read_basis(file,num);
            cache    = solver.solve(basis);
        }
        else if (action == "generate") {
            uint num    = item[1].asInt();
            real target = item[2].asDouble();
            uint trials = item[3].asInt();
            if (target < 0) {
                driver->targeting_energy = true;
                driver->target_energy    = target;
            } else {
                driver->targeting_energy = false;
                driver->target_state     = target;
            }
            driver->trial_size = trials;
            auto cd  = driver->expand_basis(basis,cache,num);
            driver->write_basis(basis, out_dir + "/basis.json");
            driver->write_convergence(cd,out_dir + "/convergence.dat",true);
        }
        else if (action == "sweep") {
            real start = item[1].asDouble();
            real end   = item[2].asDouble();
            uint steps = item[3].asInt();
            auto sd    = driver->sweep_basis(basis,start,end,steps);
            driver->write_sweep(sd,out_dir+"/sweep.dat",true);
        }
        else {
            std::cout << "Queue action not recognised.\n";
        }
    }

    /*real step_size = pi/6 / 50;
    for (real theta=0; theta<pi/6; theta += step_size) {
        SweepData sd = driver->sweep_basis(basis,theta,pi/6,1);
        driver->write_sweep(sd,out_dir+"/sweep.dat",true);
    }*/

    delete system;
    delete sample_space;
    delete driver;
}

void init(std::string file, System* system, Driver* driver, SampleSpace* sample_space)
{
    Json::Value  root;
    Json::Reader reader;

    std::ifstream config_doc(file, std::ifstream::binary);
    config_doc >> root;

    //============ PARSE SYSTEM =============
    std::cout << "Parsing System."; std::flush(std::cout);

    short particle_id_counter = 0;

    //Parse particles
    const Json::Value j_particles = root["particles"];
    for (uint n = 0; n<j_particles.size(); ++n) {
        uint count = j_particles[n].get("count",1).asInt();
        for (uint i=0; i<count; ++i) {
            ParticleType ptype = ParticleType::PT_None;
            if (j_particles[n].get("type","").asString() == "boson")
                ptype = ParticleType::PT_Boson;
            if (j_particles[n].get("type","").asString() == "fermion")
                ptype = ParticleType::PT_Fermion;

            Particle p(ptype);

            p.id           = particle_id_counter++;
            p.name         = j_particles[n].get("name","").asString();
            p.mass         = j_particles[n].get("mass",1).asDouble();
            p.identicality = j_particles[n].get("identicality",-10).asInt();
            if (p.identicality == -10) {
                std::cout << "Particle " << p.name << " has unspecified identicality";
                throw;
            }

            system->add_particle(p);
        }
    }

    std::cout << "."; std::flush(std::cout);

    //Parse interaction potentials
    const Json::Value j_interactions = root["interactions"];
    for (uint k=0; k<j_interactions.size(); ++k) {
        std::string p1 = j_interactions[k]["pair"][0].asString();
        std::string p2 = j_interactions[k]["pair"][1].asString();
        std::string t  = j_interactions[k].get("type","none").asString();

        if (t == "gaussian") {
            Interaction V;
            V.type = Interaction::Type::Gaussian;
            V.v0   = j_interactions[k].get("V0",1).asDouble();
            V.r0sq = j_interactions[k].get("r0",1).asDouble();
            V.r0sq = V.r0sq * V.r0sq;
            system->set_interaction(p1,p2,V);
        }
        else if (t == "multigaussian") {
            Interaction V;
            V.type = Interaction::Type::MultiGaussian;

            const Json::Value j_multV0 = j_interactions[k]["V0"];
            const Json::Value j_multr0 = j_interactions[k]["r0"];

            assert (j_multV0.size() == j_multr0.size());

            for (uint l=0; l<j_multV0.size(); ++l) {
                V.mult_v0.push_back( j_multV0[l].asDouble() );
                V.mult_r0sq.push_back( j_multr0[l].asDouble()
                                     * j_multr0[l].asDouble() );
            }

            system->set_interaction(p1,p2,V);
        }
        else if (t == "powerlaw") {
            Interaction V;
            V.type = Interaction::Type::PowerLaw;
            V.v0   = j_interactions[k].get("V0",1).asDouble();
            V.pow = j_interactions[k].get("pow",1).asDouble();
            system->set_interaction(p1,p2,V);
        }
        else if (t == "multipower") {
            Interaction V;
            V.type = Interaction::Type::MultiPower;

            const Json::Value j_multV0  = j_interactions[k]["V0"];
            const Json::Value j_multpow = j_interactions[k]["pow"];

            assert (j_multV0.size() == j_multpow.size());

            for (uint l=0; l<j_multV0.size(); ++l) {
                V.mult_v0.push_back(  j_multV0[l].asDouble() );
                V.mult_pow.push_back( j_multpow[l].asDouble()
                                    * j_multpow[l].asDouble() );
            }

            system->set_interaction(p1,p2,V);
        }
    }

    std::cout << "."; std::flush(std::cout);

    system->init();

    std::cout << ".\n";

    //========== PARSE SAMPLE SPACE =======
    std::cout << "Parsing Sample Space.\n";

    const Json::Value space   = root["sampleSpace"];
    const Json::Value dists   = space["distributions"];
    const Json::Value strains = space["strains"];

    for (uint k=0; k<strains.size(); ++k) {
        uint freq = strains[k].get("frequency",10).asInt();

        CG_Strain cgs;
        auto particles = system->get_particles();

        const Json::Value pairs = strains[k]["pairs"];
        for (uint l=0; l<pairs.size(); ++l) {
            std::string p1 = pairs[l]["pair"][0].asString();
            std::string p2 = pairs[l]["pair"][1].asString();
            std::string d  = pairs[l].get("distribution","").asString();

            short id1 = -1;
            short id2 = -1;

            for (auto p : particles) {
                if (p.name == p1) id1 = p.id;
                if (p.name == p2) id2 = p.id;
            }

            if (id1 == -1) std::cout << "ERROR: No particle with name" << p1 << "\n";
            if (id2 == -1) std::cout << "ERROR: No particle with name" << p2 << "\n";

            for (uint i=0; i<dists.size(); ++i) {
                if (dists[i].get("name","none").asString() == d) {
                    if (dists[i].get("type","").asString() == "uniform") {
                        double min = dists[i].get("min",0).asDouble();
                        double max = dists[i].get("max",0).asDouble();

                        cgs.set_distribution(id1,id2,
                                             std::shared_ptr<SamplingDistribution>(
                                                 new SD_Uniform(min,max)));
                    }
                    else if (dists[i].get("type","").asString() == "gaussian") {
                        std::string name = dists[i].get("name","").asString();

                        double mean  = dists[i].get("mean",0).asDouble();
                        double std   = dists[i].get("std",0).asDouble();
                        double mstdf = dists[i].get("mstdf",0).asDouble();
                        double min   = dists[i].get("min",-1).asDouble();
                        double max   = dists[i].get("max",-1).asDouble();
                        double learn = dists[i].get("learn",0).asInt();
                        double hist  = dists[i].get("history",0).asInt();

                        bool has_min = false;
                        bool has_max = false;
                        if (min != -1) has_min = true;
                        if (max != -1) has_max = true;

                        cgs.set_distribution(id1,id2,
                                             std::shared_ptr<SamplingDistribution>(
                                                 new SD_Gaussian(name,mean,std,mstdf,has_min,
                                                                 min,has_max,max,learn,
                                                                 hist)));
                    }
                    else {
                        std::cout << "ERROR: distribution type invalid for "
                                  << dists[i].get("name","") << "\n";
                        throw;
                    }
                }
            }
        }

        sample_space->add_strain(cgs,freq);
    }

    //========== PARSE DRIVER =============
    std::cout << "Parsing Driver.\n";

    const Json::Value JDriver  = root["driver"];
    driver->target_state       = JDriver.get("targetState",0).asInt();
    //driver->targetEnergy     = JDriver.get("targetEnergy",1111).asDouble();
    driver->trial_size         = JDriver.get("trialSize",1).asInt();
    driver->threads            = JDriver.get("threads",1).asInt();
    driver->singularity_limit  = JDriver.get("singularityLimit",1e-10).asDouble();
    //driver->forceDiversity   = JDriver.get("forceDiversity",0).asInt();
}
