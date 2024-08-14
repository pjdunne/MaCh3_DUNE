#include <TROOT.h>
#include "samplePDFDUNEAtmBase.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"

samplePDFDUNEAtmBase::samplePDFDUNEAtmBase(double pot, std::string mc_version, covarianceXsec* xsec_cov) : samplePDFFDBase(pot, mc_version, xsec_cov) { 
}

samplePDFDUNEAtmBase::~samplePDFDUNEAtmBase() {
}

void samplePDFDUNEAtmBase::Init() {
  //Bools
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  //Inputs
  std::string mtupleprefix = SampleManager->raw()["InputFiles"]["mtupleprefix"].as<std::string>();
  std::string mtuplesuffix = SampleManager->raw()["InputFiles"]["mtuplesuffix"].as<std::string>();
  std::string splineprefix = SampleManager->raw()["InputFiles"]["splineprefix"].as<std::string>();
  std::string splinesuffix = SampleManager->raw()["InputFiles"]["splinesuffix"].as<std::string>();

  std::vector<std::string> mtuple_files;
  std::vector<std::string> spline_files;
  std::vector<int> sample_vecno;
  std::vector<int> sample_oscnutype;
  std::vector<int> sample_nutype;
  std::vector<bool> sample_signal;
  
  //Loop over all the sub-samples
  for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
    std::cout << "Found sub sample" << std::endl;
    mtuple_files.push_back(osc_channel["mtuplefile"].as<std::string>());
    spline_files.push_back(osc_channel["splinefile"].as<std::string>());
    sample_vecno.push_back(osc_channel["samplevecno"].as<int>());
    sample_nutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["nutype"].as<int>())));
    sample_oscnutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["oscnutype"].as<int>())));
    sample_signal.push_back(osc_channel["signal"].as<bool>());
  }
  
  //Now loop over the kinematic cuts
  for ( auto const &SelectionCuts : SampleManager->raw()["SelectionCuts"]) {
    std::cout << "Looping over selection cuts " << std::endl;
    SelectionStr.push_back(SelectionCuts["KinematicStr"].as<std::string>());
    
    SelectionBounds.push_back(SelectionCuts["Bounds"].as<std::vector<double>>());
    std::cout << "Found cut on string " << SelectionCuts["KinematicStr"].as<std::string>() << std::endl;
    std::cout << "With bounds " << SelectionCuts["Bounds"].as<std::vector<double>>()[0] << " to " << SelectionCuts["Bounds"].as<std::vector<double>>()[1] << std::endl;
  }
  NSelections = SelectionStr.size();
  
  // create dunemc storage
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }
  
  //Now down with yaml file for sample
  delete SampleManager;
  std::cout << "Oscnutype size: " << sample_oscnutype.size() << ", dunemcSamples size: " << dunemcSamples.size() << std::endl;  
  if(sample_oscnutype.size() != dunemcSamples.size()){std::cerr << "[ERROR:] samplePDFDUNEAtmBase::samplePDFDUNEAtmBase() - something went wrong either getting information from sample config" << std::endl; throw;}
  
  for(unsigned iSample=0 ; iSample < dunemcSamples.size() ; iSample++){
    setupDUNEMC((mtupleprefix+mtuple_files[iSample]+mtuplesuffix).c_str(), &dunemcSamples[sample_vecno[iSample]], pot, sample_nutype[iSample], sample_oscnutype[iSample], sample_signal[iSample]);
  }
  
  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    setupFDMC(&dunemcSamples[sample_vecno[iSample]], &MCSamples[sample_vecno[iSample]]);
  }
  
  std::cout << "################" << std::endl;
  std::cout << "Setup FD MC   " << std::endl;
  std::cout << "################" << std::endl;

  std::vector<std::string> spline_filepaths;

  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    spline_filepaths.push_back(splineprefix+spline_files[iSample]+splinesuffix);
  }

  splineFile = new splinesDUNE(XsecCov);

  //////////////////////////////////
  // Now add samples to spline monolith
  //////////////////////////////////

  //ETA - do we need to do this here?
  //Can't we have one splineFile object for all samples as Dan does in atmospherics fit?
  //Then just add spline files to monolith?
  std::cout<<"Adding samples to spline monolith"<<std::endl;
  std::cout << "samplename is " << samplename << std::endl;
  std::cout << "BinningOpt is " << BinningOpt << std::endl;
  std::cout << "SampleDetID is " << SampleDetID << std::endl;
  std::cout << "spline_filepaths is of size " << spline_filepaths.size() << std::endl;

  splineFile->AddSample(samplename, BinningOpt, SampleDetID, spline_filepaths);

  // Print statements for debugging
  splineFile->PrintArrayDimension();
  splineFile->CountNumberOfLoadedSplines(false, 1);
  splineFile->TransferToMonolith();
  std::cout << "--------------------------------" <<std::endl;

  std::cout << "################" << std::endl;
  std::cout << "Setup FD splines   " << std::endl;
  std::cout << "################" << std::endl;

  SetupNormParameters();
  SetupWeightPointers();

  fillSplineBins();

  _sampleFile->Close();

  std::cout << "-------------------------------------------------------------------" <<std::endl;
}


void samplePDFDUNEAtmBase::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
	for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
	  //DB Setting total weight pointers
	  MCSamples[i].ntotal_weight_pointers[j] = 6;
	  MCSamples[i].total_weight_pointers[j] = new double*[MCSamples[i].ntotal_weight_pointers[j]];
	  MCSamples[i].total_weight_pointers[j][0] = &(dunemcSamples[i].pot_s);
	  MCSamples[i].total_weight_pointers[j][1] = &(dunemcSamples[i].norm_s);
	  MCSamples[i].total_weight_pointers[j][2] = &(MCSamples[i].osc_w[j]);
	  MCSamples[i].total_weight_pointers[j][3] = &(dunemcSamples[i].rw_berpaacvwgt[j]);
	  MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].flux_w[j]);
	  MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
	}
  }

}


void samplePDFDUNEAtmBase::setupDUNEMC(const char *sampleFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats) {
  // set up splines
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;
  
  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("cafTree");
  // _data = (TTree*)_sampleFile->Get("caf");

  if(_data){
    std::cout << "Found mtuple tree is " << sampleFile << std::endl;
    std::cout << "N of entries: " << _data->GetEntries() << std::endl;
  }
  
  _data->SetBranchStatus("*", 0);
  _data->SetBranchStatus("Ev", 1);
  _data->SetBranchAddress("Ev", &_ev);
  _data->SetBranchStatus("Ev_reco_numu", 1);
  _data->SetBranchAddress("Ev_reco_numu", &_erec);
  _data->SetBranchStatus("Ev_reco_nue", 1);
  _data->SetBranchAddress("Ev_reco_nue", &_erec_nue);
  _data->SetBranchStatus("mode",1);
  _data->SetBranchAddress("mode",&_mode);
  _data->SetBranchAddress("NuMomX", &_NuMomX);
  _data->SetBranchStatus("NuMomY", 1);
  _data->SetBranchAddress("NuMomY", &_NuMomY);
  _data->SetBranchAddress("NuMomZ", &_NuMomZ);
  _data->SetBranchStatus("LepMomX", 1);
  _data->SetBranchAddress("LepMomX", &_LepMomX);
  _data->SetBranchStatus("LepMomY", 1);
  _data->SetBranchAddress("LepMomY", &_LepMomY);
  _data->SetBranchStatus("LepMomZ", 1);
  _data->SetBranchAddress("LepMomZ", &_LepMomZ);
  _data->SetBranchStatus("cvnnumu",1);
  _data->SetBranchAddress("cvnnumu", &_cvnnumu);
  _data->SetBranchStatus("cvnnue",1);
  _data->SetBranchAddress("cvnnue", &_cvnnue);
  _data->SetBranchStatus("isCC", 1);
  _data->SetBranchAddress("isCC", &_isCC);
  _data->SetBranchStatus("nuPDGunosc", 1);
  _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
  _data->SetBranchStatus("nuPDG", 1);
  _data->SetBranchAddress("nuPDG", &_nuPDGunosc);
  _data->SetBranchStatus("run", 1);
  _data->SetBranchAddress("run", &_run);
  _data->SetBranchStatus("isFHC", 1);
  _data->SetBranchAddress("isFHC", &_isFHC);
  _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
  _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);
  _data->SetBranchStatus("vtx_x", 1);
  _data->SetBranchAddress("vtx_x", &_vtx_x);
  _data->SetBranchStatus("vtx_y", 1);
  _data->SetBranchAddress("vtx_y", &_vtx_y);
  _data->SetBranchStatus("vtx_z", 1);
  _data->SetBranchAddress("vtx_z", &_vtx_z);  
  _data->SetBranchStatus("LepPDG",1); 
  _data->SetBranchAddress("LepPDG",&_ipnu); 
  _data->SetBranchStatus("LepNuAngle",1); 
  _data->SetBranchAddress("LepNuAngle",&_LepTheta); 
  _data->SetBranchStatus("Q2",1); 
  _data->SetBranchAddress("Q2",&_Q2);
  _data->SetBranchStatus("weight",1); 
  _data->SetBranchAddress("weight",&_weight);

  duneobj->norm_s = 1;
  duneobj->pot_s = 1;

  std::cout<< "pot_s = " << duneobj->pot_s << std::endl;
  std::cout<< "norm_s = " << duneobj->norm_s << std::endl;

  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;
  
  std::cout << "signal: " << duneobj->signal << std::endl;
  std::cout << "nevents: " << duneobj->nEvents << std::endl;

  // allocate memory for dunemc variables
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->rw_Q2 = new double[duneobj->nEvents];
  duneobj->beam_w = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu = new double[duneobj->nEvents];
  duneobj->rw_cvnnue = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_run = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_truecz = new double[duneobj->nEvents];

  duneobj->mode = new int[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];

  _data->GetEntry(0);
  
  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) // Loop through tree
    {
      _data->GetEntry(i);

      if (iselike) {
        duneobj->rw_erec[i] = _erec_nue;
      } else {
	duneobj->rw_erec[i] = _erec;
      }

      duneobj->rw_etru[i] = _ev; // in GeV

      duneobj->rw_cvnnumu[i] = _cvnnumu; 
      duneobj->rw_cvnnue[i] = _cvnnue;
      duneobj->rw_theta[i] = _LepMomY/sqrt(pow(_LepMomZ, 2) + pow(_LepMomY, 2) + pow(_LepMomX, 2));
      duneobj->rw_Q2[i] = _Q2;
      duneobj->rw_isCC[i] = _isCC;
      duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
      duneobj->rw_nuPDG[i] = _nuPDG;
      duneobj->rw_run[i] = _run;
      duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;
      duneobj->rw_vtx_x[i] = _vtx_x;
      duneobj->rw_vtx_y[i] = _vtx_y;
      duneobj->rw_vtx_z[i] = _vtx_z;
      duneobj->rw_truecz[i] = _NuMomY/_ev;

      //Assume everything is on Argon for now....
      duneobj->Target[i] = 40;
 
      duneobj->beam_w[i] = 1.0;
      int mode= TMath::Abs(_mode);       

      duneobj->mode[i]=SIMBMode_ToMaCh3Mode(mode, _isCC);
      duneobj->flux_w[i] = _weight;
    }
  
  _sampleFile->Close();
  std::cout << "Sample set up OK" << std::endl;
}

double samplePDFDUNEAtmBase::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
 KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter)); 
 double KinematicValue = -999;

 switch(KinPar){
 case kTrueNeutrinoEnergy:
   KinematicValue = dunemcSamples[iSample].rw_etru[iEvent]; 
   break;
 case kTrueXPos:
   KinematicValue = dunemcSamples[iSample].rw_vtx_x[iEvent];
   break;
 case kTrueYPos:
   KinematicValue = dunemcSamples[iSample].rw_vtx_y[iEvent];
   break;
 case kTrueZPos:
   KinematicValue = dunemcSamples[iSample].rw_vtx_z[iEvent];
   break;
 case kTrueCosZ:
   KinematicValue = dunemcSamples[iSample].rw_truecz[iEvent];
   break;
 default:
   std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type '" << KinematicParameter << "'" << std::endl;
   throw;
 }
 
 return KinematicValue;
}

double samplePDFDUNEAtmBase::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
  double KinematicValue = -999;

  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = dunemcSamples[iSample].rw_etru[iEvent]; 
    break;
  case kTrueXPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_x[iEvent];
    break;
  case kTrueYPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_y[iEvent];
    break;
  case kTrueZPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_z[iEvent];
    break;
  case kTrueCosZ:
    KinematicValue = dunemcSamples[iSample].rw_truecz[iEvent];
    break;
  default:
    std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
    throw;
  }

  return KinematicValue;
}

void samplePDFDUNEAtmBase::setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj)  {
  fdobj->nEvents = duneobj->nEvents;
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->x_var = new double*[fdobj->nEvents];
  fdobj->y_var = new double*[fdobj->nEvents];
  fdobj->rw_etru = new double*[fdobj->nEvents];
  fdobj->XBin = new int[fdobj->nEvents];
  fdobj->YBin = new int[fdobj->nEvents];    
  fdobj->NomXBin = new int[fdobj->nEvents];
  fdobj->NomYBin = new int[fdobj->nEvents];
  fdobj->XBin = new int [fdobj->nEvents];
  fdobj->YBin = new int [fdobj->nEvents];;	 
  fdobj->rw_lower_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_lower_lower_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_upper_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_upper_upper_xbinedge = new double [fdobj->nEvents];
  fdobj->mode = new int*[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents]; 
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->nxsec_norm_pointers = new int[fdobj->nEvents];
  fdobj->xsec_norm_pointers = new const double**[fdobj->nEvents];
  fdobj->xsec_norms_bins = new std::list< int >[fdobj->nEvents];
  fdobj->xsec_w = new double[fdobj->nEvents];
  fdobj->flux_w = new double[fdobj->nEvents];
  fdobj->osc_w = new double[fdobj->nEvents];
  fdobj->isNC = new bool[fdobj->nEvents];
  fdobj->rw_etru = new double*[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents];
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->ntotal_weight_pointers = new int[fdobj->nEvents];
  fdobj->total_weight_pointers = new double**[fdobj->nEvents];
  fdobj->Target = new int*[fdobj->nEvents];

  //DB Example Atmospheric Implementation
  fdobj->Unity = 1.;
  fdobj->osc_w_pointer = new const double*[fdobj->nEvents];
  fdobj->rw_truecz = new double[fdobj->nEvents];

  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
	fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
	fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
	fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
	//ETA - set these to a dummy value to begin with, these get set in set1DBinning or set2DBinning
	fdobj->NomXBin[iEvent] = -1;
	fdobj->NomYBin[iEvent] = -1;
	fdobj->XBin[iEvent] = -1;
	fdobj->YBin[iEvent] = -1;	 
	fdobj->rw_lower_xbinedge[iEvent] = -1;
	fdobj->rw_lower_lower_xbinedge[iEvent] = -1;
	fdobj->rw_upper_xbinedge[iEvent] = -1;
	fdobj->rw_upper_upper_xbinedge[iEvent] = -1;
	fdobj->xsec_w[iEvent] = 1.0;
	fdobj->osc_w[iEvent] = 1.0;
	fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
	fdobj->flux_w[iEvent] = duneobj->flux_w[iEvent];
	fdobj->SampleDetID = SampleDetID;

	//DB Example Atmospheric Implementation
	fdobj->osc_w_pointer[iEvent] = &(fdobj->Unity);
	fdobj->rw_truecz[iEvent] = duneobj->rw_truecz[iEvent];

	//ETA - this is where the variables that you want to bin your samples in are defined
	//If you want to bin in different variables this is where you put it for now
	switch(BinningOpt){
	case 0:
	case 1:
	  //Just point to xvar to the address of the variable you want to bin in
	  //This way we don't have to update both fdmc and skmc when we apply shifts
	  //to variables we're binning in
	  fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
	  fdobj->y_var[iEvent] = &(duneobj->rw_truecz[iEvent]);
	  break;
	case 2:
	  //Just point to xvar to the address of the variable you want to bin in
	  //This way we don't have to update both fdmc and skmc when we apply shifts
	  //to variables we're binning in
	  fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
	  fdobj->y_var[iEvent] = &(duneobj->rw_theta[iEvent]);
	  break;
	default:
	  std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__ << " unrecognised binning option" << BinningOpt << std::endl;
	  throw;
	  break;
	}
	
  }
}

std::vector<double> samplePDFDUNEAtmBase::ReturnKinematicParameterBinning(std::string KinematicParameterStr)  {
  std::cout << "ReturnKinematicVarBinning" << std::endl;
  std::vector<double> binningVector;
  return binningVector;
}
