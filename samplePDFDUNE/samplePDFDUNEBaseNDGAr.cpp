#include <TROOT.h>

#include "samplePDFDUNEBaseNDGAr.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"
#include "TLeaf.h"
//#include "duneanaobj/StandardRecord/"
//#define DEBUG

// Constructors for erec-binned errors

//!!rewrite execs to give arguments in new order
samplePDFDUNEBaseNDGAr::samplePDFDUNEBaseNDGAr(double pot, std::string mc_version, covarianceXsec* xsec_cov)
  : samplePDFBase(pot)
{ 
  std::cout << "- Using DUNE sample config in this file " << mc_version << std::endl;
  //ETA - safety feature so you can't pass a NULL xsec_cov
  if(xsec_cov == NULL){std::cerr << "[ERROR:] You've passed me a NULL xsec covariance matrix... I need this to setup splines!" << std::endl; throw;}
  init(pot, mc_version, xsec_cov);          

}

samplePDFDUNEBaseNDGAr::~samplePDFDUNEBaseNDGAr()
{
  delete sr;
  delete _sampleFile;
}

void samplePDFDUNEBaseNDGAr::init(double pot, std::string samplecfgfile, covarianceXsec *xsec_cov)
{

  Beta=1;
  useBeta=false;
  applyBetaNue=false;
  applyBetaDiag=false;

  //doubled_angle =true ;
  useNonDoubledAngles(true);
  if (doubled_angle) std::cout << "- Using non doubled angles for oscillation parameters" << std::endl;

  osc_binned = false;
  if (osc_binned) std::cout << "- Using binned oscillation weights" << std::endl;

  modes = new TH1D("modes","",120,-60,60);

  std::string mtupleprefix;
  std::string mtuplesuffix;
  std::string splineprefix;
  std::string splinesuffix;

  char* sample_char = (char*)samplecfgfile.c_str();
  //ETA - trying to read stuff from yaml file
  manager* SampleManager = new manager(sample_char);

  //Bools
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();


  iscalo_reco = SampleManager->raw()["SampleBools"]["iscalo_reco"].as<bool>(); //NK determine what reco used
  muonscore_threshold = SampleManager->raw()["SampleCuts"]["muonscore_threshold"].as<float>(); //NK determine what muon score threshold to use

  //Inputs
  mtupleprefix = SampleManager->raw()["InputFiles"]["mtupleprefix"].as<std::string>();
  mtuplesuffix = SampleManager->raw()["InputFiles"]["mtuplesuffix"].as<std::string>();
  splineprefix = SampleManager->raw()["InputFiles"]["splineprefix"].as<std::string>();
  splinesuffix = SampleManager->raw()["InputFiles"]["splinesuffix"].as<std::string>();

  //Binning
  BinningOpt = SampleManager->raw()["Binning"]["BinningOpt"].as<int>();  
  std::vector<double> sample_erec_bins = SampleManager->raw()["Binning"]["XVarBins"].as<std::vector<double>>();
  std::vector<double> sample_theta_bins = SampleManager->raw()["Binning"]["YVarBins"].as<std::vector<double>>();

  samplename = SampleManager->raw()["SampleName"].as<std::string>();

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

  std::vector<double> tempselection(3);
  for(int iSelec=0; iSelec<NSelections; iSelec++){
    tempselection[0] = ReturnKinematicParameterFromString(SelectionStr[iSelec]);
    tempselection[1] = SelectionBounds[iSelec][0];
    tempselection[2] = SelectionBounds[iSelec][1];
    StoredSelection.push_back(tempselection);
  }

  //Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  double erec_bin_edges[sample_erec_bins.size()];
  double theta_bin_edges[sample_theta_bins.size()];
  for(unsigned erec_i = 0 ; erec_i < sample_erec_bins.size() ; erec_i++){erec_bin_edges[erec_i] = sample_erec_bins[erec_i];}
  for(unsigned theta_i = 0 ; theta_i < sample_theta_bins.size() ; theta_i++){theta_bin_edges[theta_i] = sample_theta_bins[theta_i];}

  // create dunendgarmc storage
  int nSamples = SampleManager->raw()["NSubSamples"].as<int>();
  for (int i=0;i<nSamples;i++) {
    struct dunendgarmc_base obj = dunendgarmc_base();
    dunendgarmcSamples.push_back(obj);
  }
  //Now down with yaml file for sample
  delete SampleManager;
  std::cout << "Oscnutype size: " << sample_oscnutype.size() << ", dunendgarmcSamples size: " << dunendgarmcSamples.size() << endl;  
//  if(sample_oscnutype.size() != dunendgarmcSamples.size()){std::cerr << "[ERROR:] samplePDFDUNEBaseNDGAr::samplePDFDUNEBaseNDGAr() - something went wrong either getting information from sample config" << std::endl; throw;}

  for(unsigned iSample=0 ; iSample < dunendgarmcSamples.size() ; iSample++){
    setupDUNEMC((mtupleprefix+mtuple_files[iSample]+mtuplesuffix).c_str(), &dunendgarmcSamples[sample_vecno[iSample]], pot, sample_nutype[iSample], sample_oscnutype[iSample], sample_signal[iSample]);
  }

  for (int i=0;i<nSamples;i++) {
    struct fdmc_base obj = fdmc_base();
    MCSamples.push_back(obj);
  }

  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    setupFDMC(&dunendgarmcSamples[sample_vecno[iSample]], &MCSamples[sample_vecno[iSample]], (splineprefix+spline_files[iSample]+splinesuffix).c_str());
  }

  std::cout << "################" << std::endl;
  std::cout << "Setup FD MC   " << std::endl;
  std::cout << "################" << std::endl;

  // ETA - If xsec_cov hasn't been passed to the samplePDFDUNEBaseND constructor then it's NULL
  // and the old funcitonality is kept
  // this calls this function in the core code
  // this needs to come after setupFDMC as otherwise MCSamples.splinefile will be NULL
  setXsecCov(xsec_cov); 

  std::cout << "Now setting up Splines" << std::endl;
  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
	setupSplines(&MCSamples[sample_vecno[iSample]] , (splineprefix+spline_files[iSample]+splinesuffix).c_str(), MCSamples[iSample].nutype, MCSamples[iSample].signal);
  }

  std::cout << "################" << std::endl;
  std::cout << "Setup FD splines   " << std::endl;
  std::cout << "################" << std::endl;

  setupWeightPointers();

  fillSplineBins();

  #ifdef USE_PROB3
  std::cout << "- Setup Prob3++" << std::endl;
  #else
  std::cout << "- Setup CUDAProb3" << std::endl;
  #endif
  
  //_sampleFile->Close();
  char *histname = (char*)"blah";
  char *histtitle = (char*)"blahblah";

  std::cout << "-------------------------------------------------------------------" <<std::endl;

  //The binning here is arbitrary, now we get info from cfg so the 
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D = new TH1D("hErec_nue", "Reconstructed Energy", 200, 0 , 50.0);
  dathist = new TH1D("dat_nue","",200,0, 50.0); 
  _hPDF2D = new TH2D(histname,histtitle,15,0,50.0*1000,15,0,150);
  dathist2d = new TH2D("dat2d_nue","",15,0,1500,15,0,150);

  //ETA Don't forget the -1 on the size here, as it's number of bins not bin edges
  set1DBinning(sample_erec_bins.size()-1, erec_bin_edges);
  set2DBinning(sample_erec_bins.size()-1, erec_bin_edges, sample_theta_bins.size()-1, theta_bin_edges); 

  return;
}


void samplePDFDUNEBaseNDGAr::setupWeightPointers() {
  
  for (int i = 0; i < (int)dunendgarmcSamples.size(); ++i) {
	for (int j = 0; j < dunendgarmcSamples[i].nEvents; ++j) {
	  //DB Setting total weight pointers
	  MCSamples[i].ntotal_weight_pointers[j] = 6;
	  MCSamples[i].total_weight_pointers[j] = new double*[MCSamples[i].ntotal_weight_pointers[j]];
	  MCSamples[i].total_weight_pointers[j][0] = &(dunendgarmcSamples[i].pot_s);
	  MCSamples[i].total_weight_pointers[j][1] = &(dunendgarmcSamples[i].norm_s);
	  MCSamples[i].total_weight_pointers[j][2] = &(MCSamples[i].osc_w[j]);
	  MCSamples[i].total_weight_pointers[j][3] = &(dunendgarmcSamples[i].rw_berpaacvwgt[j]);
	  MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].flux_w[j]);
	  MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
	}
  }
  std::cout<<"setup weight pointers"<<std::endl;
  return;
}


void samplePDFDUNEBaseNDGAr::setupDUNEMC(const char *sampleFile, dunendgarmc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats)
{
  
  // set up splines
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;

  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("cafTree");
  if(_data){
    std::cout << "Found mtuple tree is " << sampleFile << std::endl;
    std::cout << "N of entries: " << _data->GetEntries() << std::endl;
  }

  _data->SetBranchStatus("*", 1);


  _data->SetBranchAddress("rec", &sr);

  //FIX Commenting out for now 
  /*
  TH1D* norm = (TH1D*)_sampleFile->Get("norm");
  if(!norm){
    std::cout<< "Add a norm KEY to the root file using MakeNormHists.cxx"<<std::endl;
    std::cout << "Ignoring for now" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

*/

  // now fill the actual variables
  if (!IsRHC) 
  { 
    duneobj->norm_s = 1.0;
  }
  else 
  {
    duneobj->norm_s = 1.0;
  }

  //x10 since we're using 1/10 the MC
  duneobj->pot_s = (pot)/1e21;

  //LW - eventually add norm bins to CAFs
  //duneobj->norm_s = norm->GetBinContent(1);
  //duneobj->pot_s = pot/norm->GetBinContent(2);
  std::cout<< "pot_s = " << duneobj->pot_s << std::endl;
  std::cout<< "norm_s = " << duneobj->norm_s << std::endl;
  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;
 
  
  std::cout << "signal: " << duneobj->signal << std::endl;
  std::cout << "nevents: " << duneobj->nEvents << std::endl;

  // allocate memory for dunendgarmc variables
  duneobj->rw_yrec = new double[duneobj->nEvents];
  duneobj->rw_elep_reco = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  //duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->xsec_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 

  duneobj->energyscale_w = new double[duneobj->nEvents];
  duneobj->mode = new int[duneobj->nEvents];
  duneobj->rw_lower_erec_1d = new double[duneobj->nEvents]; //lower erec bound for bin
  duneobj->rw_upper_erec_1d = new double[duneobj->nEvents]; //upper erec bound for bin
  duneobj->rw_lower_erec_2d = new double[duneobj->nEvents]; //lower erec bound for bin
  duneobj->rw_upper_erec_2d = new double[duneobj->nEvents]; //upper erec bound for bin

  duneobj->nproton = new int[duneobj->nEvents];
  duneobj->nneutron = new int[duneobj->nEvents];
  duneobj->npip = new int[duneobj->nEvents];
  duneobj->npim = new int[duneobj->nEvents];
  duneobj->npi0 = new int[duneobj->nEvents];

  duneobj->nrecomuon = new int[duneobj->nEvents];
  duneobj->ntruemuon = new int[duneobj->nEvents];
  duneobj->nmuonsratio = new double[duneobj->nEvents];
  duneobj->ntruemuonprim = new int[duneobj->nEvents];

  duneobj->nrecoparticles = new int[duneobj->nEvents];
  duneobj->in_fdv = new bool[duneobj->nEvents];
  duneobj->rw_elep_true = new double[duneobj->nEvents];

  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_rad = new double[duneobj->nEvents];

  duneobj->rw_lep_pT = new double[duneobj->nEvents];
  duneobj->rw_lep_pZ = new double[duneobj->nEvents];

  duneobj->rw_reco_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_reco_rad = new double[duneobj->nEvents];

  //These spline bins get filled in fillSplineBins
  duneobj->enu_s_bin = new unsigned int[duneobj->nEvents];
  duneobj->erec_s_bin = new unsigned int[duneobj->nEvents];
  duneobj->yrec_s_bin = new unsigned int[duneobj->nEvents];
  duneobj->flux_bin = new int[duneobj->nEvents];
  duneobj->xsec_norms_bins = new std::list< int >[duneobj->nEvents];
  // change so only points to one
  duneobj->Target = new int[duneobj->nEvents];

  int num_no_ixns =0;
  int num_no_recparticles = 0;
  int num_in_fdv = 0;
  int num_in_fdv_noreco = 0;
  int num_notin_fdv =0;
  int num_nanenergy =0;
  int num_nanparticles =0;

  //FILL DUNE STRUCT
  for (int i = 0; i < (duneobj->nEvents); ++i) // Loop through tree
    {
     _data->GetEntry(i);
     double radius = pow((pow((sr->mc.nu[0].vtx.y+150),2) + pow((sr->mc.nu[0].vtx.z-1486),2)),0.5);
     if(std::abs(sr->mc.nu[0].vtx.x)<=209.0 &&  radius<=227.02){
       //std::cout<<"this event is within the fiducial volume"<<std::endl;
       num_in_fdv++;
       duneobj->in_fdv[i] = 1;
     }
     else{
       //std::cout<<"this event is NOT within the fiducial volume"<<std::endl;
       num_notin_fdv++;
       duneobj->in_fdv[i] = 0;
     }

     if(sr->common.ixn.ngsft == 0){
       //duneobj->rw_erec[i] = (double)(0);
       float erec_total =0;
       duneobj->rw_elep_reco[i] = double(0);
       duneobj->rw_yrec[i] = (double)(0);
       num_no_ixns++;
       duneobj->nrecoparticles[i] = (int)(0);
     }
     else{
       duneobj->nrecoparticles[i] = (int)(0);
       float erec_total =0;
       float elep_reco =0;
       float muonscore = muonscore_threshold;
       int nixns = (int)(sr->common.ixn.ngsft);
       for(int i_ixn =0; i_ixn<nixns; i_ixn++){
         int nrecoparticles = (int)(sr->common.ixn.gsft[i_ixn].part.ngsft);
         duneobj->nrecoparticles[i] += (int)(sr->common.ixn.gsft[i_ixn].part.ngsft);
         int nanparticles = 0;
         if(nrecoparticles ==0){
           double radius = pow((pow((sr->mc.nu[0].vtx.y+150.),2) + pow((sr->mc.nu[0].vtx.z-1486.),2)),0.5);
           if(std::abs(sr->mc.nu[0].vtx.x)<=209.0 || radius<=227.02){
             //std::cout<<"within fd"<<std::endl;
             num_in_fdv_noreco++;
           }   
         num_no_recparticles++;}
         for(int i_part =0; i_part<nrecoparticles; i_part++){
           float erec_part = (float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].E);
           //std::cout<<"erec_part: "<<erec_part<<std::endl;
           if(std::isnan(erec_part)){nanparticles++;}
           erec_total+=erec_part;
           if((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].score.gsft_pid.muon_score>muonscore)){
             if(erec_part>elep_reco){
               elep_reco = erec_part;
               duneobj->rw_reco_vtx_x[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.x));
               duneobj->rw_reco_vtx_y[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.y));
               duneobj->rw_reco_vtx_z[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.z));
               duneobj->rw_lep_pT[i] = (double)(pow(pow((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.x), 2) + pow((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.y), 2), 0.5));
               duneobj->rw_lep_pZ[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.z);
             }
             duneobj->nrecomuon[i]++; 
           }
         }
         num_nanparticles = num_nanparticles + (nanparticles/nrecoparticles);
       } //ADD PRIMARY LEPTON ENERGY ELEP_RECO
       if(std::isnan(erec_total)){std::cout<<"nan energy"<<std::endl; num_nanenergy++; erec_total = (float)(sr->common.ixn.gsft[0].Enu.lep_calo);}
       if(iscalo_reco){duneobj->rw_erec[i]=(double)(sr->common.ixn.gsft[0].Enu.lep_calo);  /*std::cout<<"calo erec: "<<(double)(sr->common.ixn.gsft[0].Enu.lep_calo)<<std::endl;*/}
       else{duneobj->rw_erec[i]=(double)(erec_total);}
       duneobj->rw_elep_reco[i] = (double)(elep_reco);
     }

     if(duneobj->rw_erec[i] != 0){duneobj->rw_yrec[i] = (double)(((duneobj->rw_erec[i])-(duneobj->rw_elep_reco[i]))/(duneobj->rw_erec[i]));}
     else{duneobj->rw_yrec[i] = (double)(0);}
     duneobj->rw_etru[i] = (double)(sr->mc.nu[0].E); // in GeV
     duneobj->rw_isCC[i] = (int)(sr->mc.nu[0].iscc);
     duneobj->rw_nuPDGunosc[i] = sr->mc.nu[0].pdgorig;
     duneobj->rw_nuPDG[i] = sr->mc.nu[0].pdg;
     duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;

     int ntrueparticles = (int)(sr->mc.nu[0].nprim);
     for(int i_truepart =0; i_truepart<ntrueparticles; i_truepart++){
       if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 13){
         duneobj->ntruemuon[i]++;
         duneobj->ntruemuonprim[i]++;
       }
     }
     int ntruesecparticles = (int)(sr->mc.nu[0].nsec);
     for(int i_truepart =0; i_truepart<ntruesecparticles; i_truepart++){
       if(std::abs(sr->mc.nu[0].sec[i_truepart].pdg) == 13){
         duneobj->ntruemuon[i]++;
       }
     }
      
     duneobj->nproton[i] = sr->mc.nu[0].nproton;
     duneobj->nneutron[i] = sr->mc.nu[0].nneutron;
     duneobj->npip[i] = sr->mc.nu[0].npip;
     duneobj->npim[i] = sr->mc.nu[0].npim;
     duneobj->npi0[i] = sr->mc.nu[0].npi0;

     duneobj->nmuonsratio[i] = (double)(duneobj->nrecomuon[i])/(double)(duneobj->ntruemuonprim[i]);
     duneobj->rw_vtx_x[i] = (double)(sr->mc.nu[0].vtx.x);
     duneobj->rw_vtx_y[i] = (double)(sr->mc.nu[0].vtx.y);
     duneobj->rw_vtx_z[i] = (double)(sr->mc.nu[0].vtx.z);
     
     duneobj->rw_rad[i] = (double)(pow((pow((duneobj->rw_vtx_y[i]+150),2) + pow((duneobj->rw_vtx_z[i]-1486),2)),0.5)); 
     duneobj->rw_reco_rad[i] = (double)(pow(pow((duneobj->rw_reco_vtx_y[i]+150),2) + pow((duneobj->rw_reco_vtx_z[i]-1486), 2), 0.5));
     duneobj->rw_elep_true[i] = (double)(sr->mc.nu[0].prim[0].p.E);
     //Assume everything is on Argon for now....
     duneobj->Target[i] = 40;
  
     duneobj->xsec_w[i] = 1.0;

    // fill modes
    _mode = sr->mc.nu[0].mode;
    _isCC = (int)(sr->mc.nu[0].iscc);
    modes->Fill(_mode);
    //!!possible cc1pi exception might need to be 11
    int mode= TMath::Abs(_mode);       
    duneobj->mode[i]=GENIEMode_ToMaCh3Mode(mode, _isCC);
 
    duneobj->energyscale_w[i] = 1.0;
      
    duneobj->flux_w[i] = 1.0;
    }
  std::cout<<"num evts in fdv: "<<num_in_fdv<<std::endl;
  std::cout<<"num evts not in fdv: "<<num_notin_fdv<<std::endl;
  std::cout<<"num no reco particles in fdv"<< num_in_fdv_noreco<<std::endl;
  std::cout<<"num no ixns: "<<num_no_ixns<<std::endl;
  std::cout<<"num no rec particles: "<<num_no_recparticles<<std::endl;
  std::cout<<"num nan energy: "<<num_nanenergy<<std::endl;
  std::cout<<"num nan particles: "<<num_nanparticles<<std::endl;
  std::cout << "Sample set up OK" << std::endl;
  
}

double samplePDFDUNEBaseNDGAr::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent){

 KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
 double KinematicValue = -999;
 KinematicValue = samplePDFDUNEBaseNDGAr::ReturnKinematicParameter(KinPar, iSample, iEvent);
 return KinematicValue;
}

double samplePDFDUNEBaseNDGAr::ReturnKinematicParameter(KinematicTypes KinematicParameter, int iSample, int iEvent){
 
 double KinematicValue = -999;

 switch(KinematicParameter){
   case kTrueNeutrinoEnergy:
	 KinematicValue = dunendgarmcSamples[iSample].rw_etru[iEvent]; 
	 break;
   case kPionMultiplicity:
         KinematicValue = dunendgarmcSamples[iSample].npip[iEvent]+dunendgarmcSamples[iSample].npim[iEvent+dunendgarmcSamples[iSample].npi0[iEvent]];
         break;
   case kRecoNeutrinoEnergy:
         KinematicValue = dunendgarmcSamples[iSample].rw_erec[iEvent];
         break;
   case kTrueXPos:
         KinematicValue = dunendgarmcSamples[iSample].rw_vtx_x[iEvent];
         break;
   case kTrueYPos:
         KinematicValue = dunendgarmcSamples[iSample].rw_vtx_y[iEvent];
         break;
   case kTrueZPos:
         KinematicValue = dunendgarmcSamples[iSample].rw_vtx_z[iEvent];
         break;
   case kTrueRad:
         KinematicValue = dunendgarmcSamples[iSample].rw_rad[iEvent];
         break;
   case kNRecoParticles:
         KinematicValue = dunendgarmcSamples[iSample].nrecoparticles[iEvent];
         break; 
   case kInFDV:
         KinematicValue = dunendgarmcSamples[iSample].in_fdv[iEvent];
         break;
   case kTrueMinusRecoEnergyRatio:
	 KinematicValue = (dunendgarmcSamples[iSample].rw_etru[iEvent]-dunendgarmcSamples[iSample].rw_erec[iEvent])/dunendgarmcSamples[iSample].rw_etru[iEvent];
	 break;
   case kTrueMinusRecoEnergy:
	 KinematicValue = (dunendgarmcSamples[iSample].rw_etru[iEvent]-dunendgarmcSamples[iSample].rw_erec[iEvent]);
	 break;
   case kNRecoMuons:
         KinematicValue = dunendgarmcSamples[iSample].nrecomuon[iEvent];
         break;
   case kNTrueMuons:
         KinematicValue = dunendgarmcSamples[iSample].ntruemuon[iEvent];
         break;
   case kNMuonsRecoOverTruth:
         KinematicValue = dunendgarmcSamples[iSample].nmuonsratio[iEvent];
         break;
   case kRecoLepEnergy:
         KinematicValue = dunendgarmcSamples[iSample].rw_elep_reco[iEvent];
         break;
   case kTrueLepEnergy:
         KinematicValue = dunendgarmcSamples[iSample].rw_elep_true[iEvent];
         break;
   case kRecoXPos:
         KinematicValue = dunendgarmcSamples[iSample].rw_reco_vtx_x[iEvent];
         break;
   case kRecoYPos:
         KinematicValue = dunendgarmcSamples[iSample].rw_reco_vtx_y[iEvent];
         break;
   case kRecoZPos:
         KinematicValue = dunendgarmcSamples[iSample].rw_reco_vtx_z[iEvent];
         break;
   case kRecoRad:
         KinematicValue = dunendgarmcSamples[iSample].rw_reco_rad[iEvent];
         break;
   case kLepPT:
         KinematicValue = dunendgarmcSamples[iSample].rw_lep_pT[iEvent];
         break;
   case kLepPZ:
         KinematicValue = dunendgarmcSamples[iSample].rw_lep_pZ[iEvent];
         break;
   case kM3Mode:
         KinematicValue = dunendgarmcSamples[iSample].mode[iEvent];
         break;
   default:
	 std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
	 throw;
 }

  return KinematicValue;
}


void samplePDFDUNEBaseNDGAr::setupFDMC(dunendgarmc_base *duneobj, fdmc_base *fdobj, const char *splineFile) 
{

  fdobj->nEvents = duneobj->nEvents;
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->x_var = new double*[fdobj->nEvents];
  fdobj->y_var = new double*[fdobj->nEvents];
  fdobj->enu_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->xvar_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->yvar_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->rw_etru = new double*[fdobj->nEvents];
  fdobj->XBin = new int[fdobj->nEvents];
  fdobj->YBin = new int[fdobj->nEvents];    
  fdobj->NomXBin = new int[fdobj->nEvents];
  fdobj->NomYBin = new int[fdobj->nEvents];
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
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents];
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->ntotal_weight_pointers = new int[fdobj->nEvents];
  fdobj->total_weight_pointers = new double**[fdobj->nEvents];
  fdobj->Target = new int*[fdobj->nEvents];

  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
        //std::cout<<"NK: iEvent "<<iEvent<<std::endl;
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

	//ETA - this is where the variables that you want to bin your samples in are defined
	//If you want to bin in different variables this is where you put it for now
	switch(BinningOpt){
	  case 0:
	  case 1:
		//Just point to xvar to the address of the variable you want to bin in
		//This way we don't have to update both fdmc and skmc when we apply shifts
		//to variables we're binning in
		fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
		fdobj->y_var[iEvent] = &(duneobj->dummy_y);//ETA - don't think we even need this as if we have a 1D sample we never need this, just not sure I like an unitialised variable in fdmc struct? 
		break;
	  case 2:
		//Just point to xvar to the address of the variable you want to bin in
		//This way we don't have to update both fdmc and skmc when we apply shifts
		//to variables we're binning in
		fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
		fdobj->y_var[iEvent] = &(duneobj->rw_yrec[iEvent]);
		break;
	  default:
		std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__ << " unrecognised binning option" << BinningOpt << std::endl;
		throw;
		break;
	}

  }

  return;
}

void samplePDFDUNEBaseNDGAr::setupSplines(fdmc_base *fdobj, const char *splineFile, int nutype, int signal) {

  int nevents = fdobj->nEvents;
  std::cout << "##################" << std::endl;
  std::cout << "Initialising splines from file: " << (splineFile) << std::endl;
  std::cout << "##################" << std::endl;

  switch (BinningOpt){
	case 0:
	case 1:
	  fdobj->splineFile = new splinesDUNE((char*)splineFile, nutype, nevents, fdobj->SampleDetID, xsecCov);
	  if (!(nutype==1 || nutype==-1 || nutype==2 || nutype==-2)){
		std::cerr << "problem setting up splines in erec" << std::endl;
	  }
	  break;
	case 2:
	  std::cout << "Creating splineDUNEBase" << std::endl;
	  fdobj->splineFile = new splinesDUNE((char*)splineFile, nutype, nevents, (double)BinningOpt, SampleDetID, xsecCov);
	  if (!(nutype==1 || nutype==-1 || nutype==2 || nutype==-2)) {
		std::cerr << "problem setting up splines in erec" << std::endl;
	  } 
	  break;
    default:
	  break;
  }

  // ETA - Moved SetupSplineInfoArrays to be here
  fdobj->splineFile->SetupSplineInfoArray(xsecCov);
  fdobj->splineFile->SetSplineInfoArrays();

  return;
}

//This is currently here just for show. We'll implement functional parameters soon!
double samplePDFDUNEBaseNDGAr::CalcXsecWeightFunc(int iSample, int iEvent) 
{
  return 1.0;
}
std::vector<double> samplePDFDUNEBaseNDGAr::ReturnKinematicParameterBinning(std::string KinematicParameterStr) 
{
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  return ReturnKinematicParameterBinning(KinPar);
}
std::vector<double> samplePDFDUNEBaseNDGAr::ReturnKinematicParameterBinning(KinematicTypes KinPar) 
{
  std::vector<double> binningVector;
  switch(KinPar){
    case kTrueNeutrinoEnergy:
         for(double ibins =0; ibins<10*10; ibins++){
           double binval = ibins/10;
           binningVector.push_back(binval);
         }
	 break;
    case kRecoNeutrinoEnergy:
        for(double ibins =0; ibins<10*10; ibins++){
           double binval = ibins/10;
           binningVector.push_back(binval);
         } 
         break;
    case kRecoXPos:
    case kTrueXPos:
	 for(double ibins =0; ibins<259*2; ibins++){
           binningVector.push_back(ibins-259);
         }
         break;
    case kRecoYPos:
    case kTrueYPos:
	 for(double ibins =0; ibins<277*2; ibins++){
           binningVector.push_back(ibins-277-150);
         }
	 break;
    case kRecoZPos:
    case kTrueZPos:
 	 for(double ibins =0; ibins<277*2; ibins++){
           binningVector.push_back(ibins-277+1486);
         }
	 break;
    case kPionMultiplicity:
         for(double ibins =0; ibins<10; ibins++){
           binningVector.push_back(ibins);
         }
         break;
   case kNRecoParticles:
         for(double ibins =0; ibins<50; ibins++){
           binningVector.push_back(ibins);
         }
         break; 
   case kInFDV:
         for(double ibins =0; ibins<3; ibins++){
           binningVector.push_back(ibins);
         }
         break;
   case kNMuonsRecoOverTruth:
   case kTrueMinusRecoEnergyRatio:
         for(double ibins =0; ibins<20*10; ibins++){
           binningVector.push_back(-10+(double)(ibins)/10);
         }
         break;
   case kTrueMinusRecoEnergy:
         for(double ibins =0; ibins<20*10; ibins++){
           binningVector.push_back(-10+(double)(ibins)/10);
         }
         break;
   case kNTrueMuons:
   case kNRecoMuons:
         for(double ibins =0; ibins<10; ibins++){
           binningVector.push_back(ibins);
         }
         break;
    case kRecoLepEnergy:
        for(double ibins =0; ibins<10*10; ibins++){
           binningVector.push_back((double)(ibins)/10);
         } 
         break;
    case kTrueLepEnergy:
        for(double ibins =0; ibins<10*10; ibins++){
           binningVector.push_back((double)(ibins)/10);
         } 
         break;
    case kTrueRad:
    case kRecoRad:
        for(double ibins =0; ibins<300; ibins++){
           binningVector.push_back(ibins);
        }
        break;
    case kLepPT:
    case kLepPZ:
        for(double ibins =0; ibins<10*10; ibins++){
           binningVector.push_back((double)(ibins)/10);
        }
        break;
    default:
         for(double ibins =0; ibins<10*100; ibins++){
           binningVector.push_back(ibins/100);
         }
         break;
 }

  return binningVector;
}

double samplePDFDUNEBaseNDGAr::getDiscVar(int iSample, int iEvent, int varindx) 
{
  std::cout << "getDiscVar" << std::endl;
  return 0.0;
}

double samplePDFDUNEBaseNDGAr::getCovLikelihood() 
{
  std::cout << "getCovLikelihood" << std::endl;
  return 0.0;
}

void samplePDFDUNEBaseNDGAr::printPosteriors()
{
  std::cout << "printPosteriors" << std::endl;
}


/*int samplePDFDUNEBaseNDGAr::getNMCSamples()
{
 
  return 0;
}*/

int samplePDFDUNEBaseNDGAr::getNEventsInSample(int sample)
{  
  return dunendgarmcSamples[sample].nEvents;
}

TH1D* samplePDFDUNEBaseNDGAr::get1DVarHist(std::string KinematicVar1, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill!=-1) {
    if (kChannelToFill>getNMCSamples()) {
      std::cout << "Required channel is not available. kChannelToFill should be between 0 and " << getNMCSamples() << std::endl;
      std::cout << "kChannelToFill given:" << kChannelToFill << std::endl;
      std::cout << "Exitting.." << std::endl;
      throw;
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill!=-1) {
    if (kModeToFill>kMaCh3_nModes) {
      std::cout << "Required mode is not available. kModeToFill should be between 0 and " << kMaCh3_nModes << std::endl;
      std::cout << "kModeToFill given:" << kModeToFill << std::endl;
      std::cout << "Exitting.." << std::endl;
      throw;
    }
    fMode = true;
  } else {
    fMode = false;
  }

  std::vector< std::vector<double> > SelectionVec;

  if (fMode) {
    std::vector<double> SelecMode(3);
    SelecMode[0] = kM3Mode;
    SelecMode[1] = kModeToFill;
    SelecMode[2] = kModeToFill+1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    std::vector<double> SelecChannel(3);
    SelecChannel[0] = kOscChannel;
    SelecChannel[1] = kChannelToFill;
    SelecChannel[2] = kChannelToFill+1;
    SelectionVec.push_back(SelecChannel);
  }

  return get1DVarHist(KinematicVar1,SelectionVec,WeightStyle,Axis);
}

/*! DB New version of get1DVarHist which only fills histogram with events passing IsEventSelected
 * This works by having the Selection vector, where each component of Selection is a 2 or 3 length vector
 * If Selection[i].size()==3, Selection[i][0] is the ND280KinematicType which is being cut, and only events with ND280KinematicType values between Selection[i][1] and Selection[i][2] are accepted
 */
TH1D* samplePDFDUNEBaseNDGAr::get1DVarHist(std::string KinematicVar1,std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis) {

  Selection = SelectionVec;

  for(int ibounds = 0; ibounds <SelectionVec.size(); ibounds++){
    SelectionBounds.push_back(std::vector<double>(SelectionVec[ibounds].begin()+1, SelectionVec[ibounds].end()));
    SelectionStr.push_back(ReturnKinematicParameterStringFromEnum((KinematicTypes)(SelectionVec[ibounds][0])));
  }

  for (unsigned int iStoredSelection=0;iStoredSelection<StoredSelection.size();iStoredSelection++) {
    Selection.push_back(StoredSelection[iStoredSelection]);
  }
  
/*  std::vector<std::string> SelectionStr;
  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {
    if (Selection[iSelection].size()!=3) {
      std::cerr << "Selection Vector[" << iSelection << "] is not formed correctly. Expect size == 3, given:" << Selection[iSelection].size() << std::endl;
      throw;
    }
    SelectionStr.push_back(ReturnKinematicParameterStringFromEnum((KinematicTypes)(Selection[iSelection][0])));
  }
*/
//  for(int istr =0; istr<SelectionStr.size(); istr++){std::cout<<"SelectionStr: "<<SelectionStr[istr]<<std::endl;}
//  for(int ibound =0; ibound<SelectionBounds.size(); ibound++){std::cout<<"SelectionBounds: "<<SelectionBounds[ibound][0]<<" "<<SelectionBounds[ibound][1]<<std::endl;}

  //DB Cut on OscChannel in this function due to speed increase from considering skmcSamples structure (ie. Array of length NChannels)
  bool fChannel = false;
  int kChannelToFill = -1;
  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {
    if (Selection[iSelection][0] == kOscChannel) {
      fChannel = true;
      kChannelToFill = Selection[iSelection][1];
    }
  }

  if (kChannelToFill>getNMCSamples()) {
    std::cout << "Required channel is not available. kChannelToFill should be between 0 and " << getNMCSamples() << std::endl;
    std::cout << "kChannelToFill given:" << kChannelToFill << std::endl;
    std::cout << "Exitting.." << std::endl;
    throw;
  }

  TH1D* _h1DVar;
  int kPDFBinning = kNKinematicParams; //XVar
  std::string modename = "all";
  if(SelectionVec.size() != 0){modename = std::to_string((int)(Selection[0][1]));}
  std::string HistName = KinematicVar1 +"_NDGAr_"+modename;
  if (ReturnKinematicParameterFromString(KinematicVar1)==kPDFBinning) {
    _h1DVar = (TH1D*)_hPDF1D->Clone(HistName.c_str());
    _h1DVar->Reset();
  } else {
    if (Axis) {
      _h1DVar = new TH1D(HistName.c_str(), KinematicVar1.c_str(),Axis->GetNbins(),Axis->GetXbins()->GetArray());
    } else {
      std::vector<double> xBinEdges = ReturnKinematicParameterBinning(KinematicVar1);
      _h1DVar = new TH1D(HistName.c_str(), KinematicVar1.c_str(), xBinEdges.size()-1, xBinEdges.data());
    }
  }
  _h1DVar->GetXaxis()->SetTitle(KinematicVar1.c_str());
  //This should be the same as FillArray in core basically, except that
  //events will end up in different bins
  for (int i=0;i<getNMCSamples();i++) {
    if (fChannel && (i!=kChannelToFill)) {
      continue;
    }
    for(int j=0;j<getNEventsInSample(i);j++) {

      //DB Determine which events pass selection
      if (!IsEventSelected(SelectionStr,i,j)) {
		continue;
      }
      if(MCSamples[i].isNC[j]) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
        MCSamples[i].osc_w[j] = 1.;
      }

      double Weight = getEventWeight(i,j);
	  if (WeightStyle==1) {
            Weight = 1.0;
            for (int iParam=0; iParam<MCSamples[i].ntotal_weight_pointers[j] ; ++iParam) {
	      Weight *= *(MCSamples[i].total_weight_pointers[j][iParam]);
            }
	  }
	  //ETA - not sure about this
	  //if (MCSamples[i].xsec_w[j] == 0.) continue;
          //NK changing depending on weight style
          if (WeightStyle != 0 && MCSamples[i].xsec_w[j] == 0.) continue;
          if (WeightStyle == 0){Weight = 1.0;}

	  double Var1_Val;
          KinematicTypes Var1 = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicVar1));
	  if (Var1==kPDFBinning) {
		Var1_Val = *(MCSamples[i].x_var[j]);
	  } else {
		Var1_Val = ReturnKinematicParameter(KinematicVar1,i,j);
	  }
	  if (Var1_Val!=__DEFAULT_RETURN_VAL__) {
		_h1DVar->Fill(Var1_Val,Weight);
	  }
    }
  }
  for(int iselec = 0; iselec<SelectionVec.size(); iselec++){
    SelectionBounds.pop_back();
    SelectionStr.pop_back();
  }
  /* DB: This is commented out be default
  // This code shifts the histogram meaning to Events/Bin Width but this affects the overall integral of the histogram so it should not be used anywhere we care about event rates
  // We could use Hist->Integral("width") but it would require a lot of modification throughout the code

  if (Var1!=kPDFBinning) {
    //_h1DVar->SetBinContent(1,_h1DVar->GetBinContent(0)+_h1DVar->GetBinContent(1));
    //_h1DVar->SetBinContent(_h1DVar->GetNbinsX(),_h1DVar->GetBinContent(_h1DVar->GetNbinsX())+_h1DVar->GetBinContent(_h1DVar->GetNbinsX()+1));

    for (int x=1;x<=_h1DVar->GetNbinsX();x++) {
      _h1DVar->SetBinContent(x,_h1DVar->GetBinContent(x)/_h1DVar->GetXaxis()->GetBinWidth(x));
    }

    _h1DVar->GetYaxis()->SetTitle("Events/Bin Width");
  }
  */

//  PrettifyAxis(Var1,_h1DVar->GetXaxis(),0);

  return _h1DVar;
}

TH2D* samplePDFDUNEBaseNDGAr::get2DVarHist(std::string KinematicVar1,std::string KinematicVar2, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis, TAxis* Axis2) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill!=-1) {
    if (kChannelToFill>getNMCSamples()) {
      std::cout << "Required channel is not available. kChannelToFill should be between 0 and " << getNMCSamples() << std::endl;
      std::cout << "kChannelToFill given:" << kChannelToFill << std::endl;
      std::cout << "Exitting.." << std::endl;
      throw;
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill!=-1) {
    if (kModeToFill>kMaCh3_nModes) {
      std::cout << "Required mode is not available. kModeToFill should be between 0 and " << kMaCh3_nModes << std::endl;
      std::cout << "kModeToFill given:" << kModeToFill << std::endl;
      std::cout << "Exitting.." << std::endl;
      throw;
    }
    fMode = true;
  } else {
    fMode = false;
  }

  std::vector< std::vector<double> > SelectionVec;

  if (fMode) {
    std::vector<double> SelecMode(3);
    SelecMode[0] = kM3Mode;
    SelecMode[1] = kModeToFill;
    SelecMode[2] = kModeToFill+1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    std::vector<double> SelecChannel(3);
    SelecChannel[0] = kOscChannel;
    SelecChannel[1] = kChannelToFill;
    SelecChannel[2] = kChannelToFill+1;
    SelectionVec.push_back(SelecChannel);
  }

  return get2DVarHist(KinematicVar1,KinematicVar2, SelectionVec,WeightStyle,Axis, Axis2);
}

/*! DB New version of get1DVarHist which only fills histogram with events passing IsEventSelected
 * This works by having the Selection vector, where each component of Selection is a 2 or 3 length vector
 * If Selection[i].size()==3, Selection[i][0] is the ND280KinematicType which is being cut, and only events with ND280KinematicType values between Selection[i][1] and Selection[i][2] are accepted
 */
TH2D* samplePDFDUNEBaseNDGAr::get2DVarHist(std::string KinematicVar1,std::string KinematicVar2, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis, TAxis* Axis2) {

  Selection = SelectionVec;

  for(int ibounds = 0; ibounds <SelectionVec.size(); ibounds++){
    SelectionBounds.push_back(std::vector<double>(SelectionVec[ibounds].begin()+1, SelectionVec[ibounds].end()));
    SelectionStr.push_back(ReturnKinematicParameterStringFromEnum((KinematicTypes)(SelectionVec[ibounds][0])));
  }

  for (unsigned int iStoredSelection=0;iStoredSelection<StoredSelection.size();iStoredSelection++) {
    Selection.push_back(StoredSelection[iStoredSelection]);
  }
  
  //DB Cut on OscChannel in this function due to speed increase from considering skmcSamples structure (ie. Array of length NChannels)
  bool fChannel = false;
  int kChannelToFill = -1;
  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {
    if (Selection[iSelection][0] == kOscChannel) {
      fChannel = true;
      kChannelToFill = Selection[iSelection][1];
    }
  }

  if (kChannelToFill>getNMCSamples()) {
    std::cout << "Required channel is not available. kChannelToFill should be between 0 and " << getNMCSamples() << std::endl;
    std::cout << "kChannelToFill given:" << kChannelToFill << std::endl;
    std::cout << "Exitting.." << std::endl;
    throw;
  }

  TH2D* _h2DVar;
  int kPDFBinning = kNKinematicParams; //XVar
  std::string modename = "all";
  if(SelectionVec.size() != 0){modename = std::to_string((int)(Selection[0][1]));}
  std::string HistName = KinematicVar1 +"_"+KinematicVar2+"_NDGAr_"+modename;
  if (ReturnKinematicParameterFromString(KinematicVar1)==kPDFBinning) {
    _h2DVar = (TH2D*)_hPDF2D->Clone(HistName.c_str());
    _h2DVar->Reset();
  } else {
    if (Axis) {
      _h2DVar = new TH2D(HistName.c_str(), KinematicVar1.c_str(),Axis->GetNbins(),Axis->GetXbins()->GetArray(), Axis2->GetNbins(),Axis2->GetXbins()->GetArray());
    } else {
      std::vector<double> xBinEdges = ReturnKinematicParameterBinning(KinematicVar1);
      std::vector<double> yBinEdges = ReturnKinematicParameterBinning(KinematicVar2);
      _h2DVar = new TH2D(HistName.c_str(), KinematicVar1.c_str(), xBinEdges.size()-1, xBinEdges.data(), yBinEdges.size()-1, yBinEdges.data());
    }
  }
  _h2DVar->GetXaxis()->SetTitle(KinematicVar1.c_str());
  _h2DVar->GetYaxis()->SetTitle(KinematicVar2.c_str());
   //This should be the same as FillArray in core basically, except that
  //events will end up in different bins
  for (int i=0;i<getNMCSamples();i++) {
    if (fChannel && (i!=kChannelToFill)) {
      continue;
    }
    for(int j=0;j<getNEventsInSample(i);j++) {

      //DB Determine which events pass selection
      if (!IsEventSelected(SelectionStr,i,j)) {
		continue;
      }
      if(MCSamples[i].isNC[j]) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
        MCSamples[i].osc_w[j] = 1.;
      }

      double Weight = getEventWeight(i,j);
	  if (WeightStyle==1) {
            Weight = 1.0;
            for (int iParam=0; iParam<MCSamples[i].ntotal_weight_pointers[j] ; ++iParam) {
	      Weight *= *(MCSamples[i].total_weight_pointers[j][iParam]);
            }
 	  }
	  //ETA - not sure about this
	  //if (MCSamples[i].xsec_w[j] == 0.) continue;
          //NK changing depending on weight style
          if (WeightStyle != 0 && MCSamples[i].xsec_w[j] == 0.) continue;
          if (WeightStyle == 0){Weight = 1.0;}

	  double Var1_Val;
          double Var2_Val;
          KinematicTypes Var1 = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicVar1));
          KinematicTypes Var2 = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicVar2));
	  if (Var1==kPDFBinning) {
		Var1_Val = *(MCSamples[i].x_var[j]);
	  } else {
		Var1_Val = ReturnKinematicParameter(KinematicVar1,i,j);
	  }
	  if (Var2==kPDFBinning) {
		Var2_Val = *(MCSamples[i].y_var[j]);
	  } else {
		Var2_Val = ReturnKinematicParameter(KinematicVar2,i,j);
	  }
	  if (Var1_Val!=__DEFAULT_RETURN_VAL__ && Var2_Val!=__DEFAULT_RETURN_VAL__ ) {
		_h2DVar->Fill(Var1_Val,Var2_Val,Weight);
	  }
    }
  }
  for(int iselec = 0; iselec<SelectionVec.size(); iselec++){
    SelectionBounds.pop_back();
    SelectionStr.pop_back();
  }
  /* DB: This is commented out be default
  // This code shifts the histogram meaning to Events/Bin Width but this affects the overall integral of the histogram so it should not be used anywhere we care about event rates
  // We could use Hist->Integral("width") but it would require a lot of modification throughout the code

  if (Var1!=kPDFBinning) {
    //_h1DVar->SetBinContent(1,_h1DVar->GetBinContent(0)+_h1DVar->GetBinContent(1));
    //_h1DVar->SetBinContent(_h1DVar->GetNbinsX(),_h1DVar->GetBinContent(_h1DVar->GetNbinsX())+_h1DVar->GetBinContent(_h1DVar->GetNbinsX()+1));

    for (int x=1;x<=_h1DVar->GetNbinsX();x++) {
      _h1DVar->SetBinContent(x,_h1DVar->GetBinContent(x)/_h1DVar->GetXaxis()->GetBinWidth(x));
    }

    _h1DVar->GetYaxis()->SetTitle("Events/Bin Width");
  }
  */

//  PrettifyAxis(Var1,_h1DVar->GetXaxis(),0);

  return _h2DVar;
}
