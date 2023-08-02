#include <atmoscpupropagator.hpp>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <fenv.h>

using FLOAT_T = double;

int main(){
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    std::cout << std::setprecision(10);

    std::vector<FLOAT_T> cosineList;
    cosineList.push_back(-0.99895);

    std::vector<FLOAT_T> energyList;
    energyList.push_back(0.0035);

    int n_cosines = 1;
    int n_energies = 1;

    int nu_flav = -1;

    const FLOAT_T theta12 = 0.5695951908800630486710466089860865317151404697548723;
    const FLOAT_T theta13 = 0.1608752771983210967007023071793306595103776477788280;
    const FLOAT_T theta23 = 0.7853981633974483096156608458198757210492923498437764;
    const FLOAT_T dcp     = 0.0;

    const FLOAT_T dm12sq = 7.9e-5;
    const FLOAT_T dm23sq = 2.5e-3;

    int n_threads = 1;
    cudaprob3::AtmosCpuPropagator<FLOAT_T> *mypropagator;
    mypropagator = new cudaprob3::AtmosCpuPropagator<FLOAT_T>(n_cosines, n_energies, n_threads); // cpu propagator with 4 threads

    mypropagator->setEnergyList(energyList);
    mypropagator->setCosineList(cosineList);
    mypropagator->setDensityFromFile("/dune/app/users/barrowd/CUDAProb3/MaCh3_DUNE/build/_deps/cudaprob3-src/models/PREM_4layer.dat");
    mypropagator->setMNSMatrix(theta12, theta13, theta23, dcp, nu_flav);
    mypropagator->setNeutrinoMasses(dm12sq, dm23sq);
    mypropagator->setProductionHeight(22.0);
    mypropagator->calculateProbabilities(cudaprob3::Antineutrino);
    std::cout << mypropagator->getProbability(0,0, cudaprob3::ProbType::e_e) << std::endl;
}
