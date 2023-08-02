#include "TCanvas.h"
#include <atmoscudapropagator.cu> // include openmp propagator
#include <hpc_helpers.cuh> // timer

#include <atmoscpupropagator.hpp>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include "TH2F.h"
#include "TStyle.h"

#define PI 3.14159265359


// using namespace cudaprob3; // namespace of the propagators
using FLOAT_T = float;

struct oscParams
{
  FLOAT_T theta12 = 0.5695951908800630486710466089860865317151404697548723;
  FLOAT_T theta13 = 0.1608752771983210967007023071793306595103776477788280;
  FLOAT_T theta23 = 0.7853981633974483096156608458198757210492923498437764;
  FLOAT_T dcp     = 0.0;
  FLOAT_T dm12sq = 7.9e-5;
  FLOAT_T dm23sq = 2.5e-3;
};

struct flavour
{
  int pdg;
  std::string name;
  std::string name_latex;
};

struct oscType
{
  FLOAT_T dcp = 0.;
  int hierarchy = 1;
  cudaprob3::NeutrinoType nuType = cudaprob3::NeutrinoType::Neutrino;
};

std::map<int, struct flavour> flavours = {
  {12, {12, "nue", "#nu_{e}"}},
  {14, {14, "numu", "#nu_{#mu}"}},
  {16, {16, "nutau", "#nu_{#tau}"}},
  {-12, {-12, "nuebar", "#bar{#nu}_{e}"}},
  {-14, {-14, "numubar", "#bar{#nu}_{#mu}"}},
  {-16, {-16, "nutaubar", "#bar{#nu}_{#tau}"}}
};

std::map<cudaprob3::ProbType, std::pair<int, int>> transitions = {
  {cudaprob3::ProbType::e_e, {12, 12}},
  {cudaprob3::ProbType::e_m, {12, 14}},
  {cudaprob3::ProbType::e_t, {12, 16}},
  {cudaprob3::ProbType::m_e, {14, 12}},
  {cudaprob3::ProbType::m_m, {14, 14}},
  {cudaprob3::ProbType::m_t, {14, 16}}
};

std::vector<struct oscType> oscTypes = {
  {0, 1, cudaprob3::NeutrinoType::Neutrino},
  {0, 1, cudaprob3::NeutrinoType::Antineutrino},
  {0, -1, cudaprob3::NeutrinoType::Neutrino},
  {0, -1, cudaprob3::NeutrinoType::Antineutrino},
  {-PI/2, 1, cudaprob3::NeutrinoType::Neutrino},
  {-PI/2, 1, cudaprob3::NeutrinoType::Antineutrino},
  {-PI/2, -1, cudaprob3::NeutrinoType::Neutrino},
  {-PI/2, -1, cudaprob3::NeutrinoType::Antineutrino},
};

bool operator<(const oscType& l, const oscType& r) {
     return (l.dcp<r.dcp || (l.dcp==r.dcp && l.hierarchy<r.hierarchy) || (l.dcp==r.dcp && l.hierarchy==r.hierarchy && l.nuType < r.nuType));
}

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

std::map<cudaprob3::ProbType, FLOAT_T*> get_oscillograms(int n_energies, int n_cosines, const oscParams& pars, cudaprob3::NeutrinoType nuType){
  std::map<cudaprob3::ProbType, FLOAT_T*> oscillograms;
  int n_threads = 8;
  cudaprob3::AtmosCpuPropagator<FLOAT_T> *propagator;
  propagator = new cudaprob3::AtmosCpuPropagator<FLOAT_T>(n_cosines, n_energies, n_threads); // cpu propagator with 4 threads
  std::vector<FLOAT_T> cosineList = linspace((FLOAT_T)-0.995, (FLOAT_T)1.0, n_cosines);
  std::vector<FLOAT_T> energyList = logspace((FLOAT_T)1.e-1, (FLOAT_T)1.e2, n_energies);

  propagator->setEnergyList(energyList);
  propagator->setCosineList(cosineList);

  propagator->setDensityFromFile("/home/pgranger/atmospherics/debug/MaCh3_DUNE/build/_deps/cudaprob3-src/models/PREM_4layer.dat");
  propagator->setNeutrinoMasses(pars.dm12sq, pars.dm23sq);
  propagator->setProductionHeight(22.0);
  propagator->setMNSMatrix(pars.theta12, pars.theta13, pars.theta23, pars.dcp, -1);

  propagator->calculateProbabilities(nuType);

  //first result access after calculation triggers data transfer
  for (const auto & [transition, unused] : transitions){
    FLOAT_T* array = new FLOAT_T[n_energies*n_cosines];
    propagator->getProbabilityArr(array, transition);
    oscillograms[transition] = array;
  }

  delete propagator;

  return oscillograms;
}

TH2F* make_h2(uint n_energies, uint n_cosines, FLOAT_T *data, std::string fname, std::string title){
  std::vector<FLOAT_T> cosineList = linspace((FLOAT_T)-0.995, (FLOAT_T)1.0, n_cosines);
  std::vector<FLOAT_T> energyList = logspace((FLOAT_T)1.e-1, (FLOAT_T)1.e2, n_energies);

  TH2F* h2 = new TH2F("", "", n_energies - 1, energyList.data(), n_cosines - 1, cosineList.data());
  for(uint i = 0; i < n_energies; i++){
    for(uint j = 0; j < n_cosines; j ++){
      h2->SetBinContent(i + 1, j + 1, data[i*n_cosines + j]);
    }
  }

  TCanvas c("", "");
  c.cd();
  h2->SetTitle(title.c_str());
  h2->Draw("COLZ");
  gPad->SetLogx(kTRUE);
  c.SaveAs(fname.c_str());

  return h2;
}

std::map<struct oscType, std::map<cudaprob3::ProbType, FLOAT_T*>> get_all_oscillograms(uint n_energies, uint n_cosines){
  std::map<struct oscType, std::map<cudaprob3::ProbType, FLOAT_T*>> all_oscillograms;
  for(const oscType& ot : oscTypes){
    std::cout << "dcp: " << ot.dcp << " ; hierarchy: " << ot.hierarchy << " : " << "nuType: " << ot.nuType << std::endl;
    oscParams params;
    params.dcp = ot.dcp;
    params.dm23sq *= ot.hierarchy;

    all_oscillograms[ot] = get_oscillograms(n_energies, n_cosines, params, ot.nuType);
  }

  return all_oscillograms;
}

void subtract(FLOAT_T* a, FLOAT_T* b, FLOAT_T* result, uint len){
  for(uint i = 0; i < len; i++){
    result[i] = a[i] - b[i];
  }
}

void plot_diff(const std::map<struct oscType, std::map<cudaprob3::ProbType, FLOAT_T*>> &oscillograms,
  struct oscType ref, struct oscType cmp, uint n_energies, uint n_cosines, std::string title_post, std::string fname_post){

  const std::map<cudaprob3::ProbType, FLOAT_T*>& osc_ref = oscillograms.at(ref);
  const std::map<cudaprob3::ProbType, FLOAT_T*>& osc_cmp = oscillograms.at(cmp);

  cudaprob3::NeutrinoType nuType = ref.nuType;
  int reverseFactor = (nuType == cudaprob3::NeutrinoType::Antineutrino) ? -1 : 1;

  for(const auto &[transition, ref_arr] : osc_ref){
    FLOAT_T* cmp_arr = osc_cmp.at(transition);
    FLOAT_T* diff = new FLOAT_T[n_energies*n_cosines];
    subtract(cmp_arr, ref_arr, diff, n_energies*n_cosines);

    std::pair<int, int> fl = transitions[transition];
    struct flavour init = flavours[fl.first*reverseFactor];
    struct flavour final = flavours[fl.second*reverseFactor];

    std::string fname = init.name + "_to_" + final.name + fname_post + ".root";
    std::string title = init.name_latex + " #rightarrow " + final.name_latex + title_post + ";E [GeV];cos #theta";
    make_h2(n_energies, n_cosines, diff, fname, title);

    delete diff;
  }
  
}

int main(int argc, char** argv){
  gStyle->SetOptStat(0);
  //// Binning
  int n_cosines = 100;
  int n_energies = 100;

  if(argc > 1)
    n_cosines = std::atoi(argv[1]);
  if(argc > 2)
    n_energies = std::atoi(argv[2]);

  // Prob3++ probRoot.cc parameters in radians 
  oscParams params;

  std::map<struct oscType, std::map<cudaprob3::ProbType, FLOAT_T*>> oscillograms = get_all_oscillograms(n_energies, n_cosines);

  for(cudaprob3::NeutrinoType nuType : {cudaprob3::NeutrinoType::Neutrino, cudaprob3::NeutrinoType::Antineutrino}){
    oscType nh_dcp0 = {0, 1, nuType};
    oscType ih_dcp0 = {0, -1, nuType};

    oscType nh_dcppi = {-PI/2, 1, nuType};

    plot_diff(oscillograms, nh_dcp0, ih_dcp0, n_energies, n_cosines, " (NH - IH)", "_hswap");
    plot_diff(oscillograms, nh_dcp0, nh_dcppi, n_energies, n_cosines, " (#delta_{CP} #frac{#pi}{2} - 0)", "_dscpswap");
  }

  for(auto const& [key, array] : oscillograms){
    for(auto const& [key2, arr] : array){
      delete arr;
    }
  }

  

}
