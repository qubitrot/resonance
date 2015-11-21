#include <iostream>
#include <fstream>
#include <ctime>
#include "tclap/CmdLine.h"
#include "json/json.h"
#include "typedefs.h"
#include "system.h"
#include "solver.h"
#include "driver.h"
#include "sampling.h"

void init(std::string, System*&, Solver*&, SampleSpace*&, Driver*&);
bool flag_exists(char** a, char** b, std::string& flag);

int main(int argc, char* argv[])
{
    System*      system;
    Solver*      solver;
    SampleSpace* sampleSpace;
    Driver*      driver;

    std::string outdir;
    std::string configFile;
    std::string basisFile;

    try {
        TCLAP::CmdLine cmd("Complex Scaling Method", ' ' ,"0.1");

        TCLAP::ValueArg<std::string> config("c","config","Config file",true,"","string");
        TCLAP::ValueArg<std::string> output("o","output","Output directory",false,"default","string");
        TCLAP::ValueArg<std::string> basisf("b","basis","Basis file",false,"none","string");

        cmd.add(config);
        cmd.add(output);
        cmd.add(basisf);

        cmd.parse(argc,argv);

        configFile = config.getValue();
        basisFile  = basisf.getValue();

        init(configFile,system,solver,sampleSpace,driver);

        std::string confname;

        if (output.getValue() != "default") {
            outdir = output.getValue();
        } else {
            time_t t = time(0);
            struct tm* now = localtime(&t);

            std::string date;
            date = std::to_string(now->tm_year+1900) + std::string("-")
                 + std::to_string(now->tm_mon+1)     + std::string("-")
                 + std::to_string(now->tm_mday)      + std::string("-")
                 + std::to_string(now->tm_hour);

            int lastdot   = configFile.find_last_of(".");
            int lastslash = configFile.find_last_of("/");
            confname = configFile.substr(lastslash+1,lastdot-lastslash-1);

            outdir = "data/" + date + "_" + confname;
        }

        std::string mkdir = std::string("mkdir -p ") + outdir;
        std::system(mkdir.c_str());
        std::system(("cp " + configFile + " " + outdir + "/" + confname + ".json").c_str());

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

    if (basisFile != "none") {
        driver->readBasis(basisFile,100);
    }

    std::cout << "Let's go!\n";

    for (int i=0; i<10; ++i) {
        driver->generateBasis(1);
        driver->writeBasis(outdir+"/basis.json");
        driver->writeConvergenceData(outdir+"/convergence.dat");
        //driver->targetState++;
    }

    driver->sweepAngle(0,pi/5,50);

    std::cout << "Deleting Driver.\n";
    delete driver;
    std::cout << "Deleting Sample Space.\n";
    delete sampleSpace;
    std::cout << "Deleting Solver.\n";
    delete solver;
    std::cout << "Deleting Driver.\n";
    delete system;

    return 0;
}

void init(std::string file, System*& system, Solver*& solver, SampleSpace*& space, Driver*& driver)
{
    Json::Value  root;
    Json::Reader reader;

    std::ifstream config_doc(file, std::ifstream::binary);
    config_doc >> root;

    //============ PARSE SYSTEM =============
    system = new System();

    std::cout << "Parsing System.\n";

    //Parse particles
    const Json::Value Jparticles = root["particles"];
    for (uint n = 0; n<Jparticles.size(); ++n) {
        uint count = Jparticles[n].get("count",1).asInt();
        for (uint i=0; i<count; ++i) {
            ParticleType ptype = ParticleType::PT_None;
            if (Jparticles[n].get("type","").asString() == "boson")   ptype = ParticleType::PT_Boson;
            if (Jparticles[n].get("type","").asString() == "fermion") ptype = ParticleType::PT_Fermion;

            Particle* p = new Particle(ptype);

            p->name         = Jparticles[n].get("name","").asString();
            p->mass         = Jparticles[n].get("mass",1).asDouble();
            p->identicality = Jparticles[n].get("identicality",-10).asInt();
            if (p->identicality == -10) {
                std::cout << "Particle " << p->name << " has unspecified identicality";
                throw;
            }

            system->addParticle(p);
        }
    }

    //Parse trapping potential

    //Parse interaction potentials
    const Json::Value Jinteractions = root["interactions"];
    for (uint k=0; k<Jinteractions.size(); ++k) {
        std::string p1 = Jinteractions[k]["pair"][0].asString();
        std::string p2 = Jinteractions[k]["pair"][1].asString();
        std::string t  = Jinteractions[k].get("type","none").asString();

        if (t == "gaussian") {
            InteractionV V;
            V.type = InteractionV::Type::Gaussian;
            V.v0   = Jinteractions[k].get("V0",1).asDouble();
            V.r0   = Jinteractions[k].get("r0",1).asDouble();
            system->setInteractionPotential(p1,p2,V);
        }
    }

    system->init();

    //========== PARSE SOLVER =============
    solver = new CpuSolver(system);

    std::cout << "Parsing Solver\n";

    //========== PARSE SAMPLE SPACE =======
    space = new SampleSpace();

    std::cout << "Parsing Sample Space\n";

    const Json::Value ss      = root["sampleSpace"];
    const Json::Value dists   = ss["distributions"];
    const Json::Value strains = ss["strains"];

    std::unordered_map<std::string,std::shared_ptr<SamplingDistribution> > distributions;
    for (uint k=0; k<dists.size(); ++k) {
        if (dists[k].get("type","").asString() == "uniform") {
            std::string name = dists[k].get("name","").asString();
            double min = dists[k].get("min",0).asDouble();
            double max = dists[k].get("max",0).asDouble();

            distributions[name] = std::shared_ptr<SamplingDistribution>(new SD_Uniform(min,max,k));

        } else if (dists[k].get("type","").asString() == "gaussian") {
            std::string name = dists[k].get("name","").asString();
            double min = dists[k].get("min",0).asDouble();
            double max = dists[k].get("max",0).asDouble();
            double avg = dists[k].get("mean",0).asDouble();
            double std = dists[k].get("std",1).asDouble();

            distributions[name] = std::shared_ptr<SamplingDistribution>(new SD_Gaussian(avg,std,min,max));
        }
    }

    for (uint k=0; k<strains.size(); ++k) {
        uint freq = strains[k].get("frequency",10).asInt();

        MatrixStrain* ms = new MatrixStrain(system);

        const Json::Value pairs = strains[k]["pairs"];
        for (uint l=0; l<pairs.size(); ++l) {
            std::string p1 = pairs[l]["pair"][0].asString();
            std::string p2 = pairs[l]["pair"][1].asString();
            std::string d  = pairs[l].get("distribution","").asString();

            std::shared_ptr<SamplingDistribution> dp;
            if (distributions.find(d) != distributions.end()) {
                dp = distributions[d];
            } else {
                std::cout << "ERROR: Distribution \"" << d << "\" was not specified.\n";
                throw;
            }

            ms->setDistribution(p1,p2,dp);
        }
        space->addStrain(ms,freq);
    }

    //========== PARSE DRIVER =============
    driver = new Driver(system,solver,space);

    std::cout << "Parsing Driver\n";

    const Json::Value JDriver = root["driver"];
    driver->targetState      = JDriver.get("target",0).asInt();
    driver->trialSize        = JDriver.get("trialSize",1).asInt();
    driver->numThreads       = JDriver.get("threads",1).asInt();
    driver->singularityLimit = JDriver.get("singularityLimit",5e-15).asDouble();
    driver->forceDiversity   = JDriver.get("forceDiversity",0).asInt();
}
