/*
This file is part of CUDAProb3++.

CUDAProb3++ is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CUDAProb3++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with CUDAProb3++.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <atmoscudapropagator.cu> // include openmp propagator
#include <hpc_helpers.cuh> // timer

#include <atmoscpupropagator.hpp>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

using namespace cudaprob3; // namespace of the propagators

template<class T>
std::vector<T> linspace(T Emin,T Emax,unsigned int div){
  if(div==0)
    throw std::length_error("div == 0");

  std::vector<T> linpoints(div, 0.0);

  T step_lin = (Emax - Emin)/T(div-1);

  T EE = Emin;

  for(unsigned int i=0; i<div-1; i++, EE+=step_lin)
    linpoints[i] = EE;

  linpoints[div-1] = Emax;

  return linpoints;
}

template<class T>
std::vector<T> logspace(T Emin,T Emax,unsigned int div){
  if(div==0)
    throw std::length_error("div == 0");
  std::vector<T> logpoints(div, 0.0);

  T Emin_log,Emax_log;
  Emin_log = log(Emin);
  Emax_log = log(Emax);

  T step_log = (Emax_log - Emin_log)/T(div-1);

  logpoints[0]=Emin;
  T EE = Emin_log+step_log;
  for(unsigned int i=1; i<div-1; i++, EE+=step_log)
    logpoints[i] = exp(EE);
  logpoints[div-1]=Emax;
  return logpoints;
}


int main(int argc, char** argv){

  using FLOAT_T = double;

  //// Binning
  int n_cosines = 1;
  int n_energies = 1;

  if(argc > 1)
    n_cosines = std::atoi(argv[1]);
  if(argc > 2)
    n_energies = std::atoi(argv[2]);
  //if(argc > 3)
  //	n_threads = std::atoi(argv[3]);

  //std::vector<FLOAT_T> cosineList = linspace((FLOAT_T)-0.995, (FLOAT_T)1.0, n_cosines);
  //std::vector<FLOAT_T> energyList = logspace((FLOAT_T)1.e-3, (FLOAT_T)1.e2, n_energies);

  std::vector<FLOAT_T> cosineList;
  cosineList.push_back(-0.99895);

  std::vector<FLOAT_T> energyList;
  energyList.push_back(0.00350);

  // Prob3++ probRoot.cc parameters in radians 
  const FLOAT_T theta12 = 0.5695951908800630486710466089860865317151404697548723;
  const FLOAT_T theta13 = 0.1608752771983210967007023071793306595103776477788280;
  const FLOAT_T theta23 = 0.7853981633974483096156608458198757210492923498437764;
  const FLOAT_T dcp     = 0.0;

  const FLOAT_T dm12sq = 7.9e-5;
  const FLOAT_T dm23sq = 2.5e-3;

  //Neutrino Flavour +/- = neutrino/antineutrino, 1 = e, 2 = mu, 3 = tau
  int nu_flav = 2;

/*  // NuFit 4 Params
  const FLOAT_T theta12 = asin(sqrt(0.307));
  const FLOAT_T theta13 = asin(sqrt(0.0218));
  const FLOAT_T theta23 = asin(sqrt(0.528));
  const FLOAT_T dcp     = 0.64;
  const FLOAT_T dm12sq = 7.39e-5;
  const FLOAT_T dm23sq = 2.525e-3; */

  int n_threads = 20;
  AtmosCpuPropagator<FLOAT_T> *propagator;
  propagator = new AtmosCpuPropagator<FLOAT_T>(n_cosines, n_energies, n_threads); // cpu propagator with 4 threads

  std::cout << std::fixed << std::showpoint;
  std::cout << std::setprecision(10);

  // set energy list
  propagator->setEnergyList(energyList);

  //set cosine list
  propagator->setCosineList(cosineList);

  // set density model
  //propagator->setDensityFromFile("models/PREM_12layer.dat");
  //propagator->setDensityFromFile("../models/PREM_4layer_quad_v2.dat");
  //propagator->setDensityFromFile("../models/PREM_4layer_quad_v1.dat");
  //propagator->setDensityFromFile("../models/AK135_5layer_quad.dat");
  propagator->setDensityFromFile("/dune/app/users/barrowd/CUDAProb3/MaCh3_DUNE/build/_deps/cudaprob3-src/models/PREM_4layer.dat");

  // set neutrino mass differences. unit: eV^2
  propagator->setNeutrinoMasses(dm12sq, dm23sq);

  // set neutrino production height in kilometers above earth
  propagator->setProductionHeight(22.0);

  // set mixing matrix. angles in radians
  propagator->setMNSMatrix(theta12, theta13, theta23, dcp, -1);

  // perform calculation. parameter is either cudaprob3::Neutrino or cudaprob3::Antineutrino
  auto start = std::chrono::high_resolution_clock::now();
  propagator->calculateProbabilities(cudaprob3::Antineutrino);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration  = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Time taken for " << n_energies << " energies and " << n_cosines << " cosines = " << duration.count() << " ms " << std::endl;

  //first result access after calculation triggers data transfer
  std::cout << propagator->getProbability(0,0, ProbType::e_e) << std::endl;
  throw;

#if 1
  // write output to files

  //std::ofstream outfile00("out_e_e.txt");
  //std::ofstream outfile01("out_e_m.txt");
  //std::ofstream outfile02("out_e_t.txt");
  //std::ofstream outfile10("out_m_e.txt");
  //std::ofstream outfile11("out_m_m.txt");
  //std::ofstream outfile12("out_m_t.txt");
  //std::ofstream outfile20("out_t_e.txt");
  //std::ofstream outfile21("out_t_m.txt");
  //std::ofstream outfile22("out_t_t.txt");

  //outfile00 << n_cosines << " " << n_energies <<'\n';
  //outfile01 << n_cosines << " " << n_energies <<'\n';
  //outfile02 << n_cosines << " " << n_energies <<'\n';
  //outfile10 << n_cosines << " " << n_energies <<'\n';
  //outfile11 << n_cosines << " " << n_energies <<'\n';
  //outfile12 << n_cosines << " " << n_energies <<'\n';
  //outfile20 << n_cosines << " " << n_energies <<'\n';
  //outfile21 << n_cosines << " " << n_energies <<'\n';
  //outfile22 << n_cosines << " " << n_energies <<'\n';

  for(int i = 0; i < n_cosines; i++) {
    for(int j = 0; j < n_energies; j++) {
      
      std::cout <<  propagator->getProbability(i, j, ProbType::e_e) << std::endl;
      throw;
      // ProbType::x_y is probability of transition x -> y
      //outfile00 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::e_e) << " ";
      //outfile01 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::e_m) << " ";
      //outfile02 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::e_t) << " ";
      //outfile10 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::m_e) << " ";
      ///outfile11 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::m_m) << " ";
      //outfile12 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::m_t) << " ";
      //outfile20 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::t_e) << " ";
      //outfile21 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::t_m) << " ";
      //outfile22 << std::setprecision(20) << propagator->getProbability(i, j, ProbType::t_t) << " ";
    }

    //outfile00 << '\n';
    //outfile01 << '\n';
    //outfile02 << '\n';
    //outfile10 << '\n';
    //outfile11 << '\n';
    //outfile12 << '\n';
    //outfile20 << '\n';
    //outfile21 << '\n';
    //outfile22 << '\n';
  }
  //outfile00 << '\n';
  //outfile01 << '\n';
  //outfile02 << '\n';
  //outfile10 << '\n';
  //outfile11 << '\n';
  //outfile12 << '\n';
  //outfile20 << '\n';
  //outfile21 << '\n';
  //outfile22 << '\n';

  //outfile00.flush();
  //outfile01.flush();
  //outfile02.flush();
  //outfile10.flush();
  //outfile11.flush();
  //outfile12.flush();
  //outfile20.flush();
  //outfile21.flush();
  //outfile22.flush();

#endif


}
