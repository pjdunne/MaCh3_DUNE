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

//  std::string geantprefix;
//  std::string geantsuffix;

  char* sample_char = (char*)samplecfgfile.c_str();
  //ETA - trying to read stuff from yaml file
  manager* SampleManager = new manager(sample_char);

  //Bools
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  //NK - Bools whether to read GEANT files
  incl_geant = SampleManager->raw()["SampleBools"]["incl_geant"].as<bool>();
  
  iscalo_reco = SampleManager->raw()["SampleBools"]["iscalo_reco"].as<bool>(); //NK determine what reco used
  ecal_containment = SampleManager->raw()["SampleBools"]["ecal_containment"].as<bool>(); //NK do we count containment if its stopped in ECAL
  muonscore_threshold = SampleManager->raw()["SampleCuts"]["muonscore_threshold"].as<float>(); //NK determine what muon score threshold to use
  protondEdxscore = SampleManager->raw()["SampleCuts"]["protondEdxscore_threshold"].as<float>(); //NK determine what proton score threshold to use
  protontofscore = SampleManager->raw()["SampleCuts"]["protontofscore_threshold"].as<float>();  //NK determine what muon score threshold to use
  recovertexradiusthreshold =  SampleManager->raw()["SampleCuts"]["recovertexradius_threshold"].as<float>();  //NK determine what radius threshold to use
  pionenergy_threshold = (SampleManager->raw()["SampleCuts"]["pionenergy_threshold"].as<float>())/1000; //NK determine what muon score threshold to use

  B_field = SampleManager->raw()["SampleCuts"]["B_field"].as<float>(); //NK B field value in T
  momentum_resolution_threshold = SampleManager->raw()["SampleCuts"]["momentum_resolution_threshold"].as<float>(); //NK momentum_resolution threshold, total as a fraction of momentum
  pixel_spacing = SampleManager->raw()["SampleCuts"]["pixel_spacing"].as<float>(); //NK pixel spacing in mm to find num hits in y,z plane
  spatial_resolution = SampleManager->raw()["SampleCuts"]["spatial_resolution"].as<float>(); //NK spatial resolution in mm to find  in y,z plane
  adc_sampling_frequency = SampleManager->raw()["SampleCuts"]["adc_sampling_frequency"].as<float>(); //NK sampling frequency for ADC - needed to find timing resolution and spatial resolution in x dir in MHz
  drift_velocity = SampleManager->raw()["SampleCuts"]["drift_velocity"].as<float>(); //NK drift velocity of electrons in gas - needed to find timing resolution and spatial resolution in x dir in cm/microsecond
//  average_gain = SampleManager->raw()["SampleCuts"]["average_gain"].as<float>();

  pi0_reco_efficiency = SampleManager->raw()["SampleCuts"]["pi0_reco_efficiency"].as<float>(); //efficiency for pi0 reco in ECAL
  gamma_reco_efficiency = SampleManager->raw()["SampleCuts"]["gamma_reco_efficiency"].as<float>(); //efficiency for gamma reco in ECAL

  TPCFidLength = SampleManager->raw()["SampleCuts"]["TPCFidLength"].as<double>();
  TPCFidRadius = SampleManager->raw()["SampleCuts"]["TPCFidRadius"].as<double>();
  TPCInstrumentedLength = SampleManager->raw()["SampleCuts"]["TPCInstrumentedLength"].as<double>();
  TPCInstrumentedRadius = SampleManager->raw()["SampleCuts"]["TPCInstrumentedRadius"].as<double>();
  ECALInnerRadius = SampleManager->raw()["SampleCuts"]["ECALInnerRadius"].as<double>();
  ECALOuterRadius = SampleManager->raw()["SampleCuts"]["ECALOuterRadius"].as<double>();
  ECALEndCapStart = SampleManager->raw()["SampleCuts"]["ECALEndCapStart"].as<double>();
  ECALEndCapEnd = SampleManager->raw()["SampleCuts"]["ECALEndCapEnd"].as<double>();
//  hits_per_mm = SampleManager->raw()["SampleCuts"]["hits_per_mm"].as<float>();

  //muonscore_threshold = 0.5;
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
  // NK - adding the ability to also read anatrees
  //std::vector<std::string geant_files;

  
  //Loop over all the sub-samples
  for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
	std::cout << "Found sub sample" << std::endl;
	mtuple_files.push_back(osc_channel["mtuplefile"].as<std::string>());
	spline_files.push_back(osc_channel["splinefile"].as<std::string>());
	sample_vecno.push_back(osc_channel["samplevecno"].as<int>());
	sample_nutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["nutype"].as<int>())));
	sample_oscnutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["oscnutype"].as<int>())));
	sample_signal.push_back(osc_channel["signal"].as<bool>());

//        if(incl_geant){
//          geant_files.push_back(osc_channel["geantfile"].as<std::string>());
//        }
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
    std::cout<<"SelectionStr: "<<SelectionStr[iSelec]<<std::endl;
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


void samplePDFDUNEBaseNDGAr::makePixelGrid(float pixel_spacing_cm){ //make a square pixel grid with spacing defined in yaml file. Spacing must be input in the yaml file in mm, and then is converted to cm later in the code
  int numpixelrows = floor(TPCInstrumentedRadius*2/(pixel_spacing_cm)); //find number of pixels along y and z axis.
  float centre_yboundary, centre_zboundary;
  if(numpixelrows % 2 == 0){centre_yboundary = TPC_centre_y; centre_zboundary = TPC_centre_z;}
  if(numpixelrows % 2 == 1){centre_yboundary = TPC_centre_y-(pixel_spacing_cm/2); centre_zboundary = TPC_centre_z-(pixel_spacing_cm/2);}
  std::cout<<"num pixel rows: "<<numpixelrows<<std::endl;
  pixelymin = centre_yboundary - floor(numpixelrows/2)*pixel_spacing_cm;
  pixelymax = centre_yboundary + floor(numpixelrows/2)*pixel_spacing_cm;
  pixelzmin = centre_zboundary - floor(numpixelrows/2)*pixel_spacing_cm;
  pixelzmax = centre_zboundary + floor(numpixelrows/2)*pixel_spacing_cm;

  for(int i_pixel = 0; i_pixel<numpixelrows; i_pixel++){
    yboundarypositions.push_back(pixelymin+i_pixel*pixel_spacing_cm);
    zboundarypositions.push_back(pixelzmin+i_pixel*pixel_spacing_cm);
  }
}

double samplePDFDUNEBaseNDGAr::FindNHits(float pixel_spacing_cm, float centre_circle_y, float centre_circle_z, double rad_curvature){
  //use the pixel grid method to find number of pixels hit in a track.

  int num_vertices =0; // number of vertices hit. Counting to avoid duplicating
  int num_intersections =0;

  //equation for circle = (y-y0)^2 + (z-z0)^2 = r^2

  for(int i_intersect = 1; i_intersect<=yboundarypositions.size(); i_intersect++){ //check every boundary line to see if it crossed within the TPC instrumented region
   float quadratic_ineq_y= pow(rad_curvature*100, 2)-pow((yboundarypositions[i_intersect]-centre_circle_y), 2);
   if(quadratic_ineq_y > 0){
     if((pow((pow((centre_circle_z + pow(quadratic_ineq_y, 0.5)-TPC_centre_z), 2) + pow((yboundarypositions[i_intersect]-TPC_centre_y), 2)), 0.5) <=TPCInstrumentedRadius) && (pixelzmin<(centre_circle_z + pow(quadratic_ineq_y, 0.5))<=pixelzmax)){ //check that the z coord is also on pixel plane
       num_intersections++;
     if((double)(fmod((centre_circle_z + (float)(pow(quadratic_ineq_y, 0.5)) - pixelzmin), pixel_spacing_cm)) == 0){ // this is the case when a vertex is crossed so to avoid double counting pixels
       num_vertices++;
     }
     }
     if((pow((pow((centre_circle_z - pow(quadratic_ineq_y, 0.5)-TPC_centre_z), 2) + pow((yboundarypositions[i_intersect]-TPC_centre_y), 2)), 0.5) <TPCInstrumentedRadius) && pixelzmin<(centre_circle_z - pow(quadratic_ineq_y, 0.5))<=pixelzmax){ //check that the z coord is also on pixel plane
       num_intersections++;
     if((double)(fmod((centre_circle_z - (float)(pow(quadratic_ineq_y, 0.5)) - pixelzmin), pixel_spacing_cm)) == 0){ // this is the case when a vertex is crossed so to avoid double counting pixels
       num_vertices++;
     }
     }
   }
   else if(quadratic_ineq_y == 0){ //when the pixel boundary is a tangent to the circle z = z0
      num_intersections++;
   }
  }
  for(int i_intersect = 1; i_intersect<=zboundarypositions.size(); i_intersect++){
   float quadratic_ineq_z = pow(rad_curvature*100, 2)-pow((zboundarypositions[i_intersect]-centre_circle_z), 2);
   if(quadratic_ineq_z > 0){
     if( (pow((pow((centre_circle_y + pow(quadratic_ineq_z, 0.5)-TPC_centre_y), 2) + pow((zboundarypositions[i_intersect]-TPC_centre_z), 2)), 0.5) <=TPCInstrumentedRadius) && pixelymin<(centre_circle_y + pow(quadratic_ineq_z, 0.5))<=pixelymax){ //check that the z coord is also on pixel plane
       num_intersections++;
     }
     if((pow((pow((centre_circle_y - pow(quadratic_ineq_z, 0.5)-TPC_centre_y), 2) + pow((zboundarypositions[i_intersect]-TPC_centre_z), 2)), 0.5) <=TPCInstrumentedRadius) && pixelymin<(centre_circle_y - pow(quadratic_ineq_z, 0.5))<=pixelymax){ //check that the z coord is also on pixel plane
       num_intersections++;
     }
     // already checked all vertices for duplicates before so no need to repeat that
   }
   else if(quadratic_ineq_z == 0){ //when the pixel boundary is a tangent to the circle y = y0
      num_intersections++;
   }
  }
  double N_hits = num_intersections - num_vertices + 1; //Add one for the pixel that it starts on
  return N_hits; 
}

double samplePDFDUNEBaseNDGAr::CalcBeta(double p_mag, double& bg, double& gamma){ //calculate beta (v/c)
  bg = (double)(p_mag/pdgmass); //beta*gamma
  gamma = pow((1+bg*bg), 0.5); //gamma
  double beta = bg/gamma; //beta (velocity)
  return beta;
}

double samplePDFDUNEBaseNDGAr::CalcDeDx(double beta, double bg, double gamma){ //calc de/dx in 10bar argon
  double mer = m_e/pdgmass; //electron mass/mass of particle
  double tmax = 2.*1000*m_e* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).
  //Assume tcut = tmax
  double tcut = tmax;
  //Find density effect correction (delta). Sternheimer values set in header file
  double log_bg = std::log10(bg);
  double delta = 0;
  if( log_bg >= sternheimer_X0){
    delta = 2. * std::log(10.)*log_bg - sternheimer_Cbar;
    if(log_bg < sternheimer_X1){
      delta += sternheimer_A * std::pow(sternheimer_X1 - log_bg, sternheimer_K);
    }
  }
       
  //Calculate Stopping number, B
  double B = 0.5 * std::log(2.*1000*m_e*bg*bg*tcut / (1.e-12 * pow(excitationenergy, 0.5)))- 0.5*beta*beta*(1. + tcut / tmax) - 0.5*delta;
  if(B<1.){B=1.;} //Don't let B become negative

  //Calculate dE/dX 
  double dedx = density*K_const*18*B / (39.981 * beta*beta); //18 is atomic number and 39.981 is atomic mass in g/mol
 
  return dedx;
}

void samplePDFDUNEBaseNDGAr::IsParticleAccepted(dunendgarmc_base *duneobj, int& i_truepart, int& i, int& isnotaccepted, double& highestpT, float pixel_spacing_cm, int& tot_particles){
  for(int i_anapart =0; i_anapart<_MCPStartPX->size(); i_anapart++){
    float mom_tot = pow(pow(_MCPStartPX->at(i_anapart), 2)+pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2), 0.5);
    float start_radius =pow((pow(_MCPStartY->at(i_anapart)-TPC_centre_y, 2)+pow(_MCPStartZ->at(i_anapart)-TPC_centre_z, 2)), 0.5);
    if((_PDG->at(i_anapart) == sr->mc.nu[0].prim[i_truepart].pdg) && ((double)(mom_tot) >= 0.999*(double)(sr->mc.nu[0].prim[i_truepart].p.Mag())) && ((double)(mom_tot) <= 1.001*(double)(sr->mc.nu[0].prim[i_truepart].p.Mag()))){
      //std::cout<<"tot_particles: "<<tot_particles<<std::endl;
      duneobj->particlepdg[tot_particles]=_PDG->at(i_anapart);
      duneobj->particleenergy[tot_particles]=(double)(sr->mc.nu[0].prim[i_truepart].p.E);
      duneobj->particlemomentum[tot_particles]=(double)(mom_tot);

    if(std::abs(_PDG->at(i_anapart)) == 2212 || std::abs(_PDG->at(i_anapart)) == 211 || std::abs(_PDG->at(i_anapart)) == 13 || std::abs(_PDG->at(i_anapart)) == 11 || std::abs(_PDG->at(i_anapart)) == 321){
      double bg_charged = 0; double gamma_charged = 0;
      double beta_charged = CalcBeta(mom_tot, bg_charged, gamma_charged);
      double dedx_charged = CalcDeDx(beta_charged, bg_charged, gamma_charged);
      duneobj->particlededx[tot_particles] = dedx_charged;
//      std::cout<<"dedx particle "<<tot_particles<<": "<<duneobj->particlededx[tot_particles]<<std::endl; 
    }


      double transverse_mom = pow(pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2), 0.5);
      if(transverse_mom > highestpT){
          duneobj->highestpart_pT[i] = transverse_mom;
          duneobj->highestpart_theta_angle[i] = 90 - (180/M_PI)*atan(_MCPStartPX->at(i_anapart)/transverse_mom);
          highestpT = transverse_mom;
      }
 
 
  if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 2112 && (std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 14 && (std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 12){
    bool stopsinecal_radius = false;
    bool stopsinecal_length = false;
    float end_radius = pow((pow(_MCPEndY->at(i_anapart)-TPC_centre_y, 2)+pow(_MCPEndZ->at(i_anapart)-TPC_centre_z, 2)), 0.5);
    float end_length = _MCPEndX->at(i_anapart)-TPC_centre_x;
    if(ecal_containment && end_radius>ECALInnerRadius && end_radius<ECALOuterRadius){stopsinecal_radius = true;}
    if(ecal_containment && std::abs(end_length)>ECALEndCapStart && std::abs(end_length)<ECALEndCapEnd){stopsinecal_length = true;}
  if((std::abs(end_length)>TPCInstrumentedLength && !stopsinecal_length) || (end_radius>TPCInstrumentedRadius && !stopsinecal_radius)){
    if((std::abs(_MCPStartX->at(i_anapart))-TPC_centre_x)<=TPCFidLength && start_radius<=TPCFidRadius){
      if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 13) || (std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 211) || (std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 2212) || (std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 11) || (std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 321)){
        double length_track_x;
        if(std::abs(end_length)>TPCInstrumentedLength){
          if((end_length)>=0){ length_track_x = TPCInstrumentedLength - (_MCPStartX->at(i_anapart)-TPC_centre_x);} //in cm
          else{ length_track_x = -TPCInstrumentedLength - (_MCPStartX->at(i_anapart)-TPC_centre_x);} //in cm
        }
        else{length_track_x = _MCPEndX->at(i_anapart) - _MCPStartX->at(i_anapart);} //in cm
        double length_track_y = _MCPEndY->at(i_anapart)-_MCPStartY->at(i_anapart); //in cm
        double length_track_z = _MCPEndZ->at(i_anapart)-_MCPStartZ->at(i_anapart); //in cm
        double L_yz_old = pow((pow(length_track_z, 2)+pow(length_track_y, 2)), 0.5); //in cm
        double rad_curvature = transverse_mom/(0.3*B_field); //p = 0.3*B*r where p in GeV/c, B in T, r in m
        double theta_xT = atan(_MCPStartPX->at(i_anapart)/transverse_mom); //helix is travelling in x dir as that is mag field dir
        double pitch = std::abs(2*2*rad_curvature*tan(theta_xT)); //distance between two turns of a helix in m
        double helixlength = ((std::abs(length_track_x)/100)/pitch)*pow((pow(M_PI*2*rad_curvature, 2) + pow(pitch, 2)), 0.5); //L = height/pitch*sqrt((pi*diameter)**2 + pitch**2)
        double tan_theta = tan(theta_xT);
 
        //find centre of circular path
        //need to know if its a positive or negative charge
        bool positivecharged =0;
        float centre_circle_y;
        float centre_circle_z;
        double L_yz, L_yz_chord; //length of curved track in y-z plane
        if(_PDG->at(i_anapart) == 2212 || _PDG->at(i_anapart) == 211 || _PDG->at(i_anapart) == -13 || _PDG->at(i_anapart) == -11 || _PDG->at(i_anapart) == 321){ positivecharged = 1;}
        if(positivecharged){
           centre_circle_y = _MCPStartY->at(i_anapart) + (rad_curvature*100*_MCPStartPZ->at(i_anapart)/transverse_mom); //Note plus sign here as cross product gives F in direction of ( pz j - py k) F= q v x B
           centre_circle_z = _MCPStartZ->at(i_anapart) - (rad_curvature*100*_MCPStartPY->at(i_anapart)/transverse_mom);
        }
        else if(!positivecharged){
           centre_circle_y = _MCPStartY->at(i_anapart) - (rad_curvature*100*_MCPStartPZ->at(i_anapart)/transverse_mom); //Note minus sign here as cross product gives F in direction of ( -pz j + py k)
           centre_circle_z = _MCPStartZ->at(i_anapart) + (rad_curvature*100*_MCPStartPY->at(i_anapart)/transverse_mom);
        }
        //Find Position where track leaves TPC. Intersection of two circles.
          float m_const = (TPC_centre_z - centre_circle_z)/(TPC_centre_y-centre_circle_y); //gradient of line between two intersection points
          float a_const = (pow(TPCInstrumentedRadius, 2)-pow(rad_curvature*100, 2) - (pow(TPC_centre_y, 2)-pow(centre_circle_y, 2))-(pow(TPC_centre_z, 2)-pow(centre_circle_z, 2)))/(2*(centre_circle_y-TPC_centre_y));
          float quadraticformula_b = -(2*m_const*(a_const -TPC_centre_y)+2*TPC_centre_z);
          float quadraticformula_a = pow(m_const, 2)+1;
          float quadraticformula_c = pow((a_const - TPC_centre_y), 2) +pow(TPC_centre_z, 2) - pow(TPCInstrumentedRadius,2);
 
          double z_intersect_1, y_intersect_1, z_intersect_2, y_intersect_2, z_intersect_chosen, y_intersect_chosen, theta_1, theta_2, theta_start, theta_chosen, theta_diff_1, theta_diff_2;
          if(pow(quadraticformula_b, 2)-4*quadraticformula_a*quadraticformula_c>0){
          z_intersect_1 = (-quadraticformula_b+pow((pow(quadraticformula_b, 2)-4*quadraticformula_a*quadraticformula_c), 0.5))/(2*quadraticformula_a);
          y_intersect_1 = -m_const*z_intersect_1+a_const;
          z_intersect_2 = (-quadraticformula_b-pow((pow(quadraticformula_b, 2)-4*quadraticformula_a*quadraticformula_c), 0.5))/(2*quadraticformula_a);
          y_intersect_2 = -m_const*z_intersect_2+a_const;
          theta_1 = std::abs(atan((y_intersect_1-centre_circle_y)/(z_intersect_1-centre_circle_z)));
          theta_2 = std::abs(atan((y_intersect_2-centre_circle_y)/(z_intersect_2-centre_circle_z)));
          theta_start = std::abs(atan((_MCPStartY->at(i_anapart)-centre_circle_y)/(_MCPStartZ->at(i_anapart)-centre_circle_z)));
          if((z_intersect_2-centre_circle_z)<0 && (y_intersect_2-centre_circle_y)<0){theta_2 = M_PI+theta_2;}
          else if((y_intersect_2-centre_circle_y)<0){theta_2 = 2*M_PI-theta_2;}
          else if((z_intersect_2-centre_circle_z)<0){theta_2 = M_PI - theta_2;}
          if((z_intersect_1-centre_circle_z)<0 && (y_intersect_1-centre_circle_y)<0){theta_1 = M_PI+theta_1;}
          else if((y_intersect_1-centre_circle_y)<0){theta_1 = 2*M_PI-theta_1;}
          else if((z_intersect_1-centre_circle_z)<0){theta_2 = M_PI - theta_1;}
          if((_MCPStartZ->at(i_anapart)-centre_circle_z)<0 && (_MCPStartY->at(i_anapart)-centre_circle_y)<0){theta_start = M_PI+theta_start;}
          else if((_MCPStartY->at(i_anapart)-centre_circle_y)<0){theta_start = 2*M_PI-theta_start;}
          else if((_MCPStartZ->at(i_anapart)-centre_circle_z)<0){theta_start = M_PI - theta_start;}
          if(!positivecharged){ //Lorentz force law, if positively charged, will be travelling counter clockwise around the circle. if negative charged, will travel clockwise
            theta_diff_1 = (theta_1 < theta_start) ? (theta_start - theta_1) : (2*M_PI - (theta_1 - theta_start));
            theta_diff_2 = (theta_2 < theta_start) ? (theta_start - theta_2) : (2*M_PI - (theta_2 - theta_start));
          }
          else if(positivecharged){ //Lorentz force law, if positively charged, will be travelling counter clockwise around the circle. if negative charged, will travel clockwise
            theta_diff_1 = (theta_1 > theta_start) ? (theta_1 - theta_start) : (2*M_PI - (theta_start - theta_1));
            theta_diff_2 = (theta_2 > theta_start) ? (theta_2 - theta_start) : (2*M_PI - (theta_start - theta_2));
         }
         if(theta_diff_1<theta_diff_2 && (rad_curvature*100*theta_diff_1 > (TPCInstrumentedRadius-start_radius))){theta_chosen = theta_diff_1; y_intersect_chosen = y_intersect_1; z_intersect_chosen = z_intersect_1;}
         else{theta_chosen = theta_diff_2; y_intersect_chosen = y_intersect_2; z_intersect_chosen = z_intersect_2;}             
          L_yz = rad_curvature*100*theta_chosen;
          L_yz_chord = std::abs(2*rad_curvature*100*sin(theta_chosen/2));
 //                 std::cout<<"L_yz: "<<L_yz<<" theta_chosen: "<<theta_chosen<<std::endl;
          if(std::abs(L_yz*tan_theta) > std::abs(length_track_x)){
            double nturns = std::abs(length_track_x/100)/pitch;
            double theta_intersect = fmod(nturns, 1)*2*M_PI;
 //                   L_yz = rad_curvature*100*theta_intersect;
            L_yz_chord = std::abs(2*rad_curvature*100*sin(theta_intersect/2));
          }
          else{length_track_x = std::abs(L_yz*tan_theta);}
          }
          else{
            double nturns = (length_track_x/100)/pitch;
            double theta_intersect;
            if(nturns>=1){theta_intersect = 2*M_PI;}
            else{theta_intersect = nturns*2*M_PI;} 
 //                   double theta_intersect = fmod(nturns, 1)*2*M_PI;
 //                   L_yz = rad_curvature*100*theta_intersect;
            L_yz_chord = std::abs(2*rad_curvature*100*sin(theta_intersect/2));
            theta_chosen = theta_intersect;
        }
        double N_hits = FindNHits(pixel_spacing_cm, centre_circle_y, centre_circle_z, rad_curvature);
        
        double p_mag = sr->mc.nu[0].prim[i_truepart].p.Mag();
        double bg = 0; double gamma = 0;
        double beta = CalcBeta(p_mag, bg, gamma);
        double dedx = CalcDeDx(beta, bg, gamma);

        //Calculate energy loss over this length in GeV
        double E_loss = dedx*helixlength*pow(10, -3);

        double bg_end = 0; double gamma_end = 0;
        double KE_end = sr->mc.nu[0].prim[i_truepart].p.E - pdgmass - E_loss;
        double p_mag_end = pow(pow(KE_end+pdgmass, 2)-pow(pdgmass, 2), 0.5);
        double beta_end = CalcBeta(p_mag_end, bg_end, gamma_end);
        double avg_betapT = (1/(beta_end*p_mag_end*cos(theta_xT)) + 1/(beta*transverse_mom))*0.5;
 
        double sigmax = (drift_velocity/100)*(1/(adc_sampling_frequency));
        double sigmayz = (spatial_resolution/(1000)); //needs to be in m              
        double momres_yz = transverse_mom*(pow(720/(N_hits+4), 0.5)*(sigmayz*transverse_mom/(0.3*B_field*pow(L_yz_chord/100, 2)))*pow((1-(1/21)*pow((L_yz_chord/(rad_curvature*100)), 2)), 0.5));
        double momres_ms = transverse_mom*(0.016/(0.3*B_field*(L_yz/100)*cos(theta_xT)*beta))*pow(L_yz/X0, 0.5);        
        double momres_tottransverse = pow(pow(momres_yz, 2) + pow(momres_ms, 2), 0.5);
        double momres_yz_old = transverse_mom*(pow(720/(N_hits+4), 0.5)*(sigmayz*transverse_mom/(0.3*B_field*pow(L_yz/100, 2))));
        double sigma_theta = (pow(cos(theta_xT), 2))*(pitch/(2*M_PI*rad_curvature))*pow((pow(sigmax/(std::abs(length_track_x)/100),2) +pow(momres_tottransverse/transverse_mom, 2)), 0.5);
        double momres_frac = pow(pow((momres_tottransverse/transverse_mom), 2)+pow(sigma_theta*tan_theta, 2), 0.5);
        if(momres_frac > momentum_resolution_threshold){
 
          isnotaccepted++;
          duneobj->momres_nonaccepted[i] = momres_frac;
          duneobj->pdg_nonaccepted[i] = sr->mc.nu[0].prim[i_truepart].pdg;
          duneobj->rejectedpart_theta_angle[i] = 90 - (180/M_PI)*atan(tan_theta);
          duneobj->rejectedpart_track_theta_angle[i] = std::abs((180/M_PI)*theta_chosen);
          duneobj->rejectedpart_ratioradcurvature[i] = std::abs(L_yz_chord/(100*rad_curvature));
          duneobj->rejectedpart_ptot[i] = mom_tot;
          duneobj->rejectedpart_pT[i] = transverse_mom;
          duneobj->rejectedpart_radcurvature[i] = rad_curvature;
          duneobj->rejectedpart_sigmatheta[i] = sigma_theta*tan_theta;
          duneobj->rejectedpart_sigmamom[i] = momres_yz/transverse_mom;
          duneobj->highestpart_pT[i] = transverse_mom;
          duneobj->highestpart_theta_angle[i] = 90 - (180/M_PI)*atan(tan_theta);
          duneobj->highestpart_lengthtrackx[i] = std::abs(length_track_x/100);
          duneobj->highestpart_lengthtrackyz[i] = std::abs(L_yz/100);
          duneobj->rejectedpart_beta[i] = beta;
//          if(std::abs(_PDG->at(i_anapart)) != 11 && std::abs(_PDG->at(i_anapart)) != 13 && std::abs(_PDG->at(i_anapart)) != 2212 && std::abs(_PDG->at(i_anapart)) != 211 && std::abs(_PDG->at(i_anapart)) != 321){std::cout<<" pdgmass: "<<pdgmass<<" momres_ms: "<<momres_ms<<std::endl;}
//          break;
        }
        }
        else{
          if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg)) == 111){
            double pi0recoprob =(double)(std::rand())/RAND_MAX; 
            if(pi0recoprob>pi0_reco_efficiency){
              isnotaccepted++; 
              std::cout<<"rejected pi0"<<std::endl;
              duneobj->pdg_nonaccepted[i] = sr->mc.nu[0].prim[i_truepart].pdg;
//              break;
            }
          }
          else if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg)) == 22){
            double gammarecoprob =(double)(std::rand())/RAND_MAX;
            if(gammarecoprob>gamma_reco_efficiency){
              isnotaccepted++; 
              std::cout<<"rejected gamma"<<std::endl;
              duneobj->pdg_nonaccepted[i] = sr->mc.nu[0].prim[i_truepart].pdg;
//              break;
            }
          }  
          else{
            isnotaccepted++;
            duneobj->pdg_nonaccepted[i] = sr->mc.nu[0].prim[i_truepart].pdg;
            duneobj->rejectedpart_theta_angle[i] = 90 - (180/M_PI)*tan(atan((_MCPStartPX->at(i_anapart)/(pow((pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2)), 0.5)))));
            duneobj->rejectedpart_track_theta_angle[i] = -1;
            duneobj->rejectedpart_ratioradcurvature[i] = -1;
            duneobj->rejectedpart_ptot[i] = pow(pow(_MCPStartPX->at(i_anapart), 2)+pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2)+pow(_MCPStartPX->at(i_anapart), 2), 0.5);
            duneobj->rejectedpart_pT[i] = pow(pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2)+pow(_MCPStartPX->at(i_anapart), 2), 0.5);
            duneobj->rejectedpart_radcurvature[i] = (pow((pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2)), 0.5))/(0.3*B_field);
            duneobj->rejectedpart_sigmatheta[i] = 10000;
            duneobj->rejectedpart_sigmamom[i] = 10000;
            duneobj->highestpart_theta_angle[i] = 90 - (180/M_PI)*tan(atan((_MCPStartPX->at(i_anapart)/(pow((pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2)), 0.5)))));
            duneobj->highestpart_pT[i] = pow(pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2)+pow(_MCPStartPX->at(i_anapart), 2), 0.5);
   
//            break;
          }
        }
      }
      else{
 //               std::cout<<"position not in fdv"<<std::endl;
        isnotaccepted++; duneobj->pdg_nonaccepted[i] = sr->mc.nu[0].prim[i_truepart].pdg;
//        break;
      }
    }
    else if((std::abs(end_length)>TPCInstrumentedLength && stopsinecal_length) || (end_radius>TPCInstrumentedRadius && stopsinecal_radius)){
      //std::cout<<"here"<<"i_anapart: "<<i_anapart<<" _SimHitEnergy->size(): "<<_SimHitEnergy->size()<<std::endl;
      double energydepsum = 0;
      for(int i_ecaldep=0; i_ecaldep<_SimHitEnergy->size(); i_ecaldep++){
        if(_SimHitTrkID->at(i_ecaldep) == _MCPTrkID->at(i_anapart)){energydepsum = energydepsum + _SimHitEnergy->at(i_ecaldep);}
      }
      duneobj->ecaldepositfraction[tot_particles] = energydepsum/(sr->mc.nu[0].prim[i_truepart].p.E);
      //std::cout<<"ecaldepositfraction: "<<duneobj->ecaldepositfraction[tot_particles]<<std::endl;
    }
    }
    break;
    }
  }
}

void samplePDFDUNEBaseNDGAr::setupDUNEMC(const char *sampleFile, dunendgarmc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats)
{
  
  // set up splines
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;


  std::string geantfilename(sampleFile);
  std::string prefix = "inputs/DUNE_NDGAr_CAF_files/";
  std::string::size_type i_prefix = geantfilename.find(prefix);
  if(i_prefix != std::string::npos){
    geantfilename.erase(i_prefix, prefix.length());
  }
  std::string suffix = ".root";
  std::string::size_type i_suffix = geantfilename.find(suffix);
  if(i_suffix != std::string::npos){
    geantfilename.erase(i_suffix, suffix.length());
  }
  std::string geantfilename_final = "inputs/DUNE_NDGAr_AnaTrees/"+geantfilename + "_geant.root";

  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("cafTree");
  if(_data){
    std::cout << "Found mtuple tree is " << sampleFile << std::endl;
    std::cout << "N of entries: " << _data->GetEntries() << std::endl;
  }

  if(incl_geant){
    _sampleFile_geant = new TFile(geantfilename_final.c_str(), "READ");
    _data_geant = (TTree*)_sampleFile_geant->Get("GArAnaTree");
    _data->AddFriend(_data_geant);
  }

  _data->SetBranchStatus("*", 1);


  _data->SetBranchAddress("rec", &sr);

  if(incl_geant){
    _data->SetBranchAddress("MCPStartX", &_MCPStartX);
    _data->SetBranchAddress("MCPStartY", &_MCPStartY);
    _data->SetBranchAddress("MCPStartZ", &_MCPStartZ);
    _data->SetBranchAddress("MCPEndX", &_MCPEndX);
    _data->SetBranchAddress("MCPEndY", &_MCPEndY);
    _data->SetBranchAddress("MCPEndZ", &_MCPEndZ);
    _data->SetBranchAddress("MCPStartPX", &_MCPStartPX);
    _data->SetBranchAddress("MCPStartPY", &_MCPStartPY);
    _data->SetBranchAddress("MCPStartPZ", &_MCPStartPZ);
    _data->SetBranchAddress("MCPEndPX", &_MCPEndPX);
    _data->SetBranchAddress("MCPEndPY", &_MCPEndPY);
    _data->SetBranchAddress("MCPEndPZ", &_MCPEndPZ);
    _data->SetBranchAddress("PDG", &_PDG);
    _data->SetBranchAddress("MCPTrkID", &_MCPTrkID);
    _data->SetBranchAddress("SimHitTrkID", &_SimHitTrkID);
    _data->SetBranchAddress("SimHitEnergy", &_SimHitEnergy);
  }


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

  duneobj->rw_etrurec = new double[duneobj->nEvents];
  duneobj->rw_Q2 = new double[duneobj->nEvents];
  duneobj->rw_W = new double[duneobj->nEvents];
  duneobj->rw_Q0 = new double[duneobj->nEvents];
  duneobj->rw_Q3 = new double[duneobj->nEvents];
  duneobj->rw_etrurec_nopionthreshold = new double[duneobj->nEvents];

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

  duneobj->nrecopion = new int[duneobj->nEvents];

  duneobj->nrecomuon = new int[duneobj->nEvents];
  duneobj->ntruemuon = new int[duneobj->nEvents];
  duneobj->nmuonsratio = new double[duneobj->nEvents];
  duneobj->ntruemuonprim = new int[duneobj->nEvents];
  duneobj->isnumuCC = new bool[duneobj->nEvents];

  duneobj->nrecoparticles = new int[duneobj->nEvents];
  duneobj->in_fdv = new bool[duneobj->nEvents];
  duneobj->rw_elep_true = new double[duneobj->nEvents];
  duneobj->is_accepted = new bool[duneobj->nEvents];
  duneobj->pdg_nonaccepted = new int[duneobj->nEvents];
  duneobj->momres_nonaccepted = new double[duneobj->nEvents];


  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_rad = new double[duneobj->nEvents];

  duneobj->rw_lep_pT = new double[duneobj->nEvents];
  duneobj->rw_lep_pMag = new double[duneobj->nEvents];
  duneobj->rw_lep_pZ = new double[duneobj->nEvents];
  duneobj->rw_lep_pY = new double[duneobj->nEvents];
  duneobj->rw_lep_pX = new double[duneobj->nEvents];
  duneobj->rw_reco_lep_pT = new double[duneobj->nEvents];
  duneobj->rw_reco_lep_pZ = new double[duneobj->nEvents];
  duneobj->rw_reco_lep_pY = new double[duneobj->nEvents];
  duneobj->rw_reco_lep_pX = new double[duneobj->nEvents];
  duneobj->rw_reco_lep_pMag = new double[duneobj->nEvents];
  duneobj->rw_lep_energy = new double[duneobj->nEvents];

  duneobj->rw_reco_pi_energy = new double[duneobj->nEvents];
  duneobj->rw_pi_energy = new double[duneobj->nEvents];
  duneobj->muon_pi_reco_angle = new double[duneobj->nEvents];
  duneobj->muon_pi_angle = new double[duneobj->nEvents];
  duneobj->pi_z_reco_angle = new double[duneobj->nEvents];
  duneobj->pi_z_angle = new double[duneobj->nEvents];
  duneobj->muon_z_angle = new double[duneobj->nEvents];
  duneobj->rw_reco_pi_pZ = new double[duneobj->nEvents];
  duneobj->rw_reco_pi_pY = new double[duneobj->nEvents];
  duneobj->rw_reco_pi_pX = new double[duneobj->nEvents];
  duneobj->rw_reco_pi_pT = new double[duneobj->nEvents];
  duneobj->rw_reco_pi_pMag = new double[duneobj->nEvents];
  duneobj->rw_pi_pZ = new double[duneobj->nEvents];
  duneobj->rw_pi_pY = new double[duneobj->nEvents];
  duneobj->rw_pi_pX = new double[duneobj->nEvents];
  duneobj->rw_pi_pT = new double[duneobj->nEvents];
  duneobj->rw_pi_pMag = new double[duneobj->nEvents];
  duneobj->rw_pi_min_energy = new double[duneobj->nEvents];

  duneobj->rw_reco_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_reco_rad = new double[duneobj->nEvents];

  duneobj->rejectedpart_theta_angle = new double[duneobj->nEvents];
  duneobj->rejectedpart_ptot = new double[duneobj->nEvents];
  duneobj->rejectedpart_pT = new double[duneobj->nEvents];
  duneobj->rejectedpart_radcurvature = new double[duneobj->nEvents];
  duneobj->rejectedpart_ratioradcurvature = new double[duneobj->nEvents];
  duneobj->rejectedpart_track_theta_angle = new double[duneobj->nEvents];
  duneobj->rejectedpart_sigmamom = new double[duneobj->nEvents];
  duneobj->rejectedpart_sigmatheta = new double[duneobj->nEvents];
  duneobj->highestpart_theta_angle = new double[duneobj->nEvents];
  duneobj->highestpart_pT = new double[duneobj->nEvents];
  duneobj->highestpart_lengthtrackx = new double[duneobj->nEvents];
  duneobj->highestpart_lengthtrackyz = new double[duneobj->nEvents];
  duneobj->rejectedpart_beta = new double[duneobj->nEvents];

  duneobj->ecaldepositfraction = new double[(duneobj->nEvents)*7];
  duneobj->particleevent = new int[(duneobj->nEvents)*7];
  duneobj->particlepdg = new int[(duneobj->nEvents)*7];
  duneobj->particleenergy = new double[(duneobj->nEvents)*7];
  duneobj->particlededx = new double[(duneobj->nEvents)*7];
  duneobj->particlemomentum = new double[(duneobj->nEvents)*7];

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
  int tot_particles = 0;
  int numCC = 0;
  int numFDV = 0;
  int numFDVandCC = 0;
  // NK - Make the Pixel Grid using cm as TPC coordinates are in cm

  float pixel_spacing_cm = pixel_spacing/10;
  std::srand(std::time(NULL));
  makePixelGrid(pixel_spacing_cm); //make the square pixel grid and fill two vector with the positions of the y and z pixel boundaries
  //FILL DUNE STRUCT
  for (int i = 0; i < (duneobj->nEvents); ++i) // Loop through tree
    {
     _data->GetEntry(i);
     double radius = pow((pow((sr->mc.nu[0].vtx.y-TPC_centre_y),2) + pow((sr->mc.nu[0].vtx.z-TPC_centre_z),2)),0.5); //find radius of interaction vertex
     if(std::abs(sr->mc.nu[0].vtx.x)<=TPCFidLength &&  radius<=TPCFidRadius){
//       std::cout<<"this event is within the fiducial volume"<<std::endl;
       num_in_fdv++;
       duneobj->in_fdv[i] = 1;
     }
     else{
       //std::cout<<"this event is NOT within the fiducial volume"<<std::endl;
       num_notin_fdv++;
       duneobj->in_fdv[i] = 0;
     }
     if(sr->common.ixn.ngsft == 0){ //if there is no reconstructed interaction, fill all reco variables with 0
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
       int nixns = (int)(sr->common.ixn.ngsft);
       for(int i_ixn =0; i_ixn<nixns; i_ixn++){
         int nrecoparticles = (int)(sr->common.ixn.gsft[i_ixn].part.ngsft);
         duneobj->nrecoparticles[i] += (int)(sr->common.ixn.gsft[i_ixn].part.ngsft);
         int nanparticles = 0;
         if(nrecoparticles ==0){
           double radius = pow((pow((sr->mc.nu[0].vtx.y-TPC_centre_y),2) + pow((sr->mc.nu[0].vtx.z-TPC_centre_z),2)),0.5);
           if(std::abs(sr->mc.nu[0].vtx.x)<=TPCFidLength || radius<=TPCFidRadius){
             //std::cout<<"within fd"<<std::endl;
             num_in_fdv_noreco++;
           }   
         num_no_recparticles++;}
         for(int i_part =0; i_part<nrecoparticles; i_part++){
           float erec_part = (float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].E);
           if(std::isnan(erec_part)){nanparticles++;}
           erec_total+=erec_part;
           if((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].score.gsft_pid.muon_score)>muonscore_threshold){
             if(erec_part>elep_reco){ //pick out the primary muon as the most energetic reconstructed particle with muon score > muon threshold
               elep_reco = erec_part;
               //take the reconstructed vertex as the start of the muon track
               duneobj->rw_reco_vtx_x[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.x)); 
               duneobj->rw_reco_vtx_y[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.y));
               duneobj->rw_reco_vtx_z[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.z));
               duneobj->rw_reco_lep_pT[i] = (double)(pow(pow((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.x), 2) + pow((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.y), 2), 0.5));
               duneobj->rw_reco_lep_pZ[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.z);
               duneobj->rw_reco_lep_pY[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.y);
               duneobj->rw_reco_lep_pX[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.x);
             }
             duneobj->nrecomuon[i]++; 
           }
         }
         double pionenergy =0.0;
         for(int i_part =0; i_part<nrecoparticles; i_part++){
           if(std::abs((int)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].pdg))!=0){ //particle pdg isn't 0
           if((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].score.gsft_pid.muon_score)<=muonscore_threshold){ //check the particle isn't a muon
//             std::cout<<"not a muon"<<std::endl;
             if((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].score.gsft_pid.proton_dEdx_score)<=protondEdxscore && (float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].score.gsft_pid.proton_tof_score)<=protontofscore){ //check the particle isn't a proton
//               std::cout<<"not a proton"<<std::endl;
               double partradiusvertex = pow(pow((double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.x)-duneobj->rw_reco_vtx_x[i], 2)+pow((double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.y)-duneobj->rw_reco_vtx_y[i], 2)+pow((double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.z)-duneobj->rw_reco_vtx_z[i], 2), 0.5);
               if(partradiusvertex<recovertexradiusthreshold){ // check that the particle track starts near the reco vertex
//                 std::cout<<"near reco vertex"<<std::endl;
                 if((double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].E)>pionenergy){ //pick the most energetic pion to fill the pion variables
//                   std::cout<<"energy greater than pionenergy prev"<<std::endl;                 
                   pionenergy = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].E);
                   duneobj->rw_reco_pi_energy[i] = pionenergy;
                   duneobj->rw_reco_pi_energy[i] = (double)(pionenergy-0.13957);
                   duneobj->rw_reco_pi_pZ[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.z);
                   duneobj->rw_reco_pi_pX[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.x);
                   duneobj->rw_reco_pi_pY[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.y);
                 }
                 duneobj->nrecopion[i]++;
               }
             }
           }
         }
         }
         num_nanparticles = num_nanparticles + (nanparticles/nrecoparticles);
       } //ADD PRIMARY LEPTON ENERGY ELEP_RECO
       if(std::isnan(erec_total)){
        //std::cout<<"nan energy"<<std::endl;
        num_nanenergy++;
        erec_total = (float)(sr->common.ixn.gsft[0].Enu.lep_calo);}
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
     double muonenergy =0;
     for(int i_truepart =0; i_truepart<ntrueparticles; i_truepart++){
       if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 13){
         if((double)(sr->mc.nu[0].prim[i_truepart].p.E)>muonenergy){
//           std::cout<<"muon energy: "<<(double)(sr->mc.nu[0].prim[i_truepart].p.E)<<" muon mass: 0.10566"<<std::endl;
           muonenergy = (double)(sr->mc.nu[0].prim[i_truepart].p.E);           
           duneobj->rw_lep_energy[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.E - 0.10566);
//           duneobj->rw_lep_energy[i] = (double)(pow(pow((double)(sr->mc.nu[0].prim[i_truepart].p.E), 2) - pow(0.10566, 2), 0.5));
           duneobj->rw_lep_pZ[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.pz);
           duneobj->rw_lep_pX[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.px);
           duneobj->rw_lep_pY[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.py);
         }
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
     if((bool)(sr->mc.nu[0].iscc) && (int)(sr->mc.nu[0].pdg)==14){duneobj->isnumuCC[i]=1;}
     else{duneobj->isnumuCC[i]=0;}
     duneobj->nproton[i] = sr->mc.nu[0].nproton;
     duneobj->nneutron[i] = sr->mc.nu[0].nneutron;
     duneobj->npip[i] = sr->mc.nu[0].npip;
     duneobj->npim[i] = sr->mc.nu[0].npim;
     duneobj->npi0[i] = sr->mc.nu[0].npi0;
     
//     std::cout<<"npim: "<<duneobj->npim[i]<<" npi0: "<<duneobj->npi0[i]<<" npip: "<<duneobj->npip[i]<<std::endl;
     double pionenergymax = 0;
     double pionenergymin;
     duneobj->rw_etrurec[i] = 0.0;
     int nprimpipm = 0;
     int nprimpizero =0;
     int nprimproton = 0;
     int nprimneutron = 0;
     int nprimmuon = 0;
     int isnotaccepted = 0;

     double highestpT = 0;

//     pdgmass = 0;
     duneobj->highestpart_lengthtrackx[i] = -10000;
     duneobj->highestpart_lengthtrackyz[i] = -10000;
//     if((duneobj->npip[i]+duneobj->npim[i]) == 1){
       for(int i_truepart =0; i_truepart<ntrueparticles; i_truepart++){
         if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 211){
             pdgmass = m_chargedpi;
             nprimpipm++;
           if(nprimpipm == 1){pionenergymin = (double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass); duneobj->rw_pi_min_energy[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass);}
//           std::cout<<"one pion"<<std::endl;
           if((double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass) > pionenergy_threshold){ 
             duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E);
             duneobj->rw_etrurec_nopionthreshold[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E);
             if((double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass)>pionenergymax){
               duneobj->rw_pi_energy[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass);
//           duneobj->rw_pi_energy[i] = (double)(pow(pow((double)(sr->mc.nu[0].prim[i_truepart].p.E), 2) - pow(0.13957, 2), 0.5));
               duneobj->rw_pi_pZ[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.pz);
               duneobj->rw_pi_pX[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.px);
               duneobj->rw_pi_pY[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.py);
               pionenergymax = duneobj->rw_pi_energy[i];
             }
             if((double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass)<pionenergymin && (double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass)>0){
//               std::cout<<"energy above pionenergymin: "<<(double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass)<<std::endl;
               duneobj->rw_pi_min_energy[i] = (double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass);
               pionenergymin = duneobj->rw_pi_min_energy[i];
             }
           }
           else{
             duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E - pdgmass);
             duneobj->rw_etrurec_nopionthreshold[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E); 
           }
         }
         else if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) != 211){
           if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 111){nprimpizero++; pdgmass = m_pi0; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass); duneobj->rw_etrurec_nopionthreshold[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E - pdgmass);} //duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E - 0.13498);}
           if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 2212){nprimproton++; pdgmass = m_p; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E - pdgmass); duneobj->rw_etrurec_nopionthreshold[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass);}
           if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 2112){nprimneutron++; pdgmass = m_n; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E - pdgmass); duneobj->rw_etrurec_nopionthreshold[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E-pdgmass);}
           if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 13){nprimmuon++; pdgmass = m_mu; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E); duneobj->rw_etrurec_nopionthreshold[i] += (double)(sr->mc.nu[0].prim[i_truepart].p.E);}
 
           if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 11){pdgmass = m_e;}
 
           if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 321){pdgmass = m_chargedk;}
         }
         if(isnotaccepted >0){continue;}
         duneobj->particleevent[tot_particles]= i;
         duneobj->ecaldepositfraction[tot_particles] = -1.;
         IsParticleAccepted(duneobj, i_truepart, i, isnotaccepted, highestpT, pixel_spacing_cm, tot_particles);
//         std::cout<<"dedx particle "<<tot_particles<<": "<<duneobj->particlededx[tot_particles]<<" pdg: "<<sr->mc.nu[0].prim[i_truepart].pdg<<std::endl;
         tot_particles++;
       }
       if(isnotaccepted > 0){duneobj->is_accepted[i] = 0;}
       else{duneobj->is_accepted[i] = 1;}
//     }
/*
//     PREFSI ENERGY CALC STARTS HERE
       int ntrueprefsiparticles = (int)(sr->mc.nu[0].nprefsi);
       int nprefsipipm = 0;
       int nprefsipizero =0;
       int nprefsiproton = 0;
       int nprefsineutron = 0;
       int nprefsimuon = 0;
       std::cout<<"ntrueprefsiparticles: "<<ntrueprefsiparticles<<" ntrueprimparticles: "<<ntrueparticles<<std::endl;
       for(int i_truepart =0; i_truepart<ntrueprefsiparticles; i_truepart++){
         if(std::abs(sr->mc.nu[0].prefsi[i_truepart].pdg) == 111){nprefsipizero++; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prefsi[i_truepart].p.E);}
         if(std::abs(sr->mc.nu[0].prefsi[i_truepart].pdg) == 2212){nprefsiproton++; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prefsi[i_truepart].p.E - 0.93827);}
         if(std::abs(sr->mc.nu[0].prefsi[i_truepart].pdg) == 2112){nprefsineutron++; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prefsi[i_truepart].p.E - 0.93957);}
         if(std::abs(sr->mc.nu[0].prefsi[i_truepart].pdg) == 211){nprefsipipm++; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prefsi[i_truepart].p.E);}
         if(std::abs(sr->mc.nu[0].prefsi[i_truepart].pdg) == 13){nprefsimuon++; duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prefsi[i_truepart].p.E);}
         //duneobj->rw_etrurec[i] -= (double)(sr->mc.nu[0].prefsi[i_truepart].p.E);
         //} 
//         duneobj->rw_etrurec[i] += (double)(sr->mc.nu[0].prefsi[i_truepart].p.Mag());
        std::cout<<" etrurec["<<i_truepart<<"]: "<<duneobj->rw_etrurec[i]<<std::endl; 
        }
        duneobj->rw_etrurec[i]+=(double)(sr->mc.nu[0].prim[0].p.E); //NK THIS IS FOR ADDING PRIMARY MUON ENERGY TO PREFSI PARTS
*/ 
//     PREFSI ENERGY CALC ENDS HERE

//     duneobj->rw_etrurec[i] = duneobj->rw_etrurec[i] - (nprimpipm - nprefsipipm)*0.13957 - (nprimpizero - nprefsipizero)*0.13498 - (nprimproton - nprefsiproton)*0.93827 - (nprimneutron - nprefsineutron)*0.93957; 
//     std::cout<<"true charged pi: "<<duneobj->npip[i]+duneobj->npim[i]<<" true pi0: "<<duneobj->npi0[i]<<" prefsi charged pi: "<<nprefsipipm<<" prim charged pi: "<<nprimpipm<<" prefsi pi0: "<<nprefsipizero<<" nprimpizero: "<<nprimpizero<<" true protons: "<<duneobj->nproton[i]<<" prefsi protons: "<<nprefsiproton<<" prim proton: "<<nprimproton<<" true neutrons: "<<duneobj->nneutron[i]<<" prefsi neutron: "<<nprefsineutron<<" prim neutron: "<<nprimneutron<<std::endl;
//     std::cout<<"prefsi muon: "<<nprefsimuon<<" prim muon: "<<nprimmuon<<std::endl;
//     std::cout<<"True NuE: "<<duneobj->rw_etru[i]<<"Reco NuE: "<<duneobj->rw_etrurec[i]<<" True Minus Reco: "<<(duneobj->rw_etru[i]-duneobj->rw_etrurec[i])<<" Ratio: "<<(duneobj->rw_etru[i]-duneobj->rw_etrurec[i])/duneobj->rw_etru[i]<<std::endl;

//     if((duneobj->npim[i]+duneobj->npip[i]) >0){std::cout<<"Min Pion Energy: "<<duneobj->rw_pi_min_energy[i]<<std::endl;}
     duneobj->rw_lep_pT[i] = (double)(pow(pow(duneobj->rw_lep_pX[i], 2) + pow(duneobj->rw_lep_pY[i], 2), 0.5));
     duneobj->rw_pi_pT[i] = (double)(pow(pow(duneobj->rw_pi_pX[i], 2) + pow(duneobj->rw_pi_pY[i], 2), 0.5));
     duneobj->rw_reco_pi_pT[i] = (double)(pow(pow(duneobj->rw_reco_pi_pX[i], 2) + pow(duneobj->rw_reco_pi_pY[i], 2), 0.5));
     duneobj->rw_reco_pi_pMag[i] = (double)(pow(pow(duneobj->rw_reco_pi_pX[i], 2) + pow(duneobj->rw_reco_pi_pY[i], 2) + pow(duneobj->rw_reco_pi_pZ[i], 2), 0.5));
     duneobj->rw_reco_lep_pMag[i] = (double)(pow(pow(duneobj->rw_reco_lep_pX[i], 2) + pow(duneobj->rw_reco_lep_pY[i], 2) + pow(duneobj->rw_reco_lep_pZ[i], 2), 0.5));
     duneobj->rw_pi_pMag[i] = (double)(pow(pow(duneobj->rw_pi_pX[i], 2) + pow(duneobj->rw_pi_pY[i], 2) + pow(duneobj->rw_pi_pZ[i], 2), 0.5));
     duneobj->rw_lep_pMag[i] = (double)(pow(pow(duneobj->rw_lep_pX[i], 2) + pow(duneobj->rw_lep_pY[i], 2) + pow(duneobj->rw_lep_pZ[i], 2), 0.5));
     duneobj->muon_pi_reco_angle[i] = (180/TMath::Pi())*TMath::ACos((double)((duneobj->rw_reco_lep_pX[i]*duneobj->rw_reco_pi_pX[i] + duneobj->rw_reco_lep_pY[i]*duneobj->rw_reco_pi_pY[i] + duneobj->rw_reco_lep_pZ[i]*duneobj->rw_reco_pi_pZ[i])/(duneobj->rw_reco_pi_pMag[i]*duneobj->rw_reco_lep_pMag[i])));
     duneobj->muon_pi_angle[i] = (180/TMath::Pi())*TMath::ACos((double)((duneobj->rw_lep_pX[i]*duneobj->rw_pi_pX[i] + duneobj->rw_reco_lep_pY[i]*duneobj->rw_pi_pY[i] + duneobj->rw_lep_pZ[i]*duneobj->rw_pi_pZ[i])/(duneobj->rw_pi_pMag[i]*duneobj->rw_lep_pMag[i])));
     duneobj->pi_z_reco_angle[i] = (180/TMath::Pi())*TMath::ACos((double)(duneobj->rw_reco_pi_pZ[i]/duneobj->rw_reco_pi_pMag[i]));
     duneobj->pi_z_angle[i] =  (180/TMath::Pi())*TMath::ACos((double)(duneobj->rw_pi_pZ[i]/duneobj->rw_pi_pMag[i]));
     duneobj->muon_z_angle[i] =  (180/TMath::Pi())*TMath::ACos((double)(duneobj->rw_lep_pZ[i]/duneobj->rw_lep_pMag[i]));

     if(duneobj->muon_z_angle[i]>=180){std::cout<<"muon z angle: "<<(double)(duneobj->muon_z_angle[i])<<" pZ: "<<(double)(duneobj->rw_lep_pZ[i])<<" pMag: "<<duneobj->rw_lep_pMag[i]<<std::endl;}
//     std::cout<<"reco pi energy: "<<duneobj->rw_reco_pi_energy[i]<<" true pi energu: "<<duneobj->rw_pi_energy[i]<<std::endl;
//     std::cout<<"reco pi pZ: "<<duneobj->rw_reco_pi_pZ[i]<<" reco pi pMag: "<<duneobj->rw_reco_pi_pMag[i]<<" True pi pZ: "<<duneobj->rw_pi_pZ[i]<<" True pi pMag: "<<duneobj->rw_pi_pMag[i]<<std::endl;
//     std::cout<<"reco lep Mag: "<<duneobj->rw_lep_pMag[i]<<std::endl;
     duneobj->nmuonsratio[i] = (double)(duneobj->nrecomuon[i])/(double)(duneobj->ntruemuonprim[i]);
     duneobj->rw_vtx_x[i] = (double)(sr->mc.nu[0].vtx.x);
     duneobj->rw_vtx_y[i] = (double)(sr->mc.nu[0].vtx.y);
     duneobj->rw_vtx_z[i] = (double)(sr->mc.nu[0].vtx.z);
     
     duneobj->rw_rad[i] = (double)(pow((pow((duneobj->rw_vtx_y[i]-TPC_centre_y),2) + pow((duneobj->rw_vtx_z[i]-TPC_centre_z),2)),0.5)); 
     duneobj->rw_reco_rad[i] = (double)(pow(pow((duneobj->rw_reco_vtx_y[i]-TPC_centre_y),2) + pow((duneobj->rw_reco_vtx_z[i]-TPC_centre_z), 2), 0.5));
     duneobj->rw_elep_true[i] = (double)(sr->mc.nu[0].prim[0].p.E);

     //Assume everything is on Argon for now....
     duneobj->Target[i] = 40;
  
     duneobj->xsec_w[i] = 1.0;

     duneobj->rw_W[i] = (double)(sr->mc.nu[0].W); 
     duneobj->rw_Q2[i] = (double)(sr->mc.nu[0].Q2);
     duneobj->rw_Q0[i] = (double)(sr->mc.nu[0].q0);
     duneobj->rw_Q3[i] = (double)(sr->mc.nu[0].modq);
    // fill modes
    _mode = sr->mc.nu[0].mode;
    _isCC = (int)(sr->mc.nu[0].iscc);
//    std::cout<<"mode: "<<_mode<<" IsCC: "<<_isCC<<std::endl;
    modes->Fill(_mode);
    //!!possible cc1pi exception might need to be 11
    int mode= TMath::Abs(_mode);       
    duneobj->mode[i]=GENIEMode_ToMaCh3Mode(mode, _isCC);
 
    duneobj->energyscale_w[i] = 1.0;
      
    duneobj->flux_w[i] = 1.0;
    if(duneobj->rw_isCC[i] == 1 && duneobj->in_fdv[i] == 1){
      numFDVandCC++;    
    }
    if(duneobj->rw_isCC[i] == 1){numCC++;}
    if(duneobj->in_fdv[i] == 1){numFDV++;}
  }
  std::cout<<"total particles: "<<tot_particles<<std::endl;
  std::cout<<"num events CC:"<<numCC<<std::endl;
  std::cout<<"num evenys FDV: "<<numFDV<<std::endl;
  std::cout<<"num in CC and FDV: "<<numFDVandCC<<std::endl;
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
   case kIdealNeutrinoRecoEnergy:
	 KinematicValue = dunendgarmcSamples[iSample].rw_etrurec[iEvent]; 
	 break;
   case kPionMultiplicity:
         KinematicValue = dunendgarmcSamples[iSample].npip[iEvent]+dunendgarmcSamples[iSample].npim[iEvent]+dunendgarmcSamples[iSample].npi0[iEvent];
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
   case kTrueMinusIdealRecoEnergyRatio:
	 KinematicValue = (dunendgarmcSamples[iSample].rw_etru[iEvent]-dunendgarmcSamples[iSample].rw_etrurec[iEvent])/dunendgarmcSamples[iSample].rw_etru[iEvent];
	 break;
   case kTrueMinusIdealRecoEnergy:
	 KinematicValue = (dunendgarmcSamples[iSample].rw_etru[iEvent]-dunendgarmcSamples[iSample].rw_etrurec[iEvent]);
	 break;
   case kNRecoMuons:
         KinematicValue = dunendgarmcSamples[iSample].nrecomuon[iEvent];
         break;
   case kNTruePrimMuons:
         KinematicValue = dunendgarmcSamples[iSample].ntruemuonprim[iEvent];
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
   case kIsnumuCC:
         KinematicValue = dunendgarmcSamples[iSample].isnumuCC[iEvent];
         break;
   case kM3Mode:
         KinematicValue = dunendgarmcSamples[iSample].mode[iEvent];
         break;
   case kLepRecoPT:
         KinematicValue = dunendgarmcSamples[iSample].rw_reco_lep_pT[iEvent];
         break;
   case kLepRecoPZ:
         KinematicValue = dunendgarmcSamples[iSample].rw_reco_lep_pZ[iEvent];
         break;
   case kMuonPiAngle:
         KinematicValue = dunendgarmcSamples[iSample].muon_pi_angle[iEvent];
         break;
   case kMuonPiRecoAngle:
         KinematicValue = dunendgarmcSamples[iSample].muon_pi_reco_angle[iEvent];
         break;
   case kPiRecoEnergy:
         KinematicValue = dunendgarmcSamples[iSample].rw_reco_pi_energy[iEvent];
         break;
   case kPiTrueEnergy:
         KinematicValue = dunendgarmcSamples[iSample].rw_pi_energy[iEvent];
         break;
   case kPiZAngle:
         KinematicValue = dunendgarmcSamples[iSample].pi_z_angle[iEvent];
         break;
   case kPiZRecoAngle:
         KinematicValue = dunendgarmcSamples[iSample].pi_z_reco_angle[iEvent];
         break;
   case kNChargedPions:
         KinematicValue = dunendgarmcSamples[iSample].npip[iEvent]+dunendgarmcSamples[iSample].npim[iEvent];
         break;
   case kNRecoPions:
         KinematicValue = dunendgarmcSamples[iSample].nrecopion[iEvent];
         break;
   case kPiRecoMomentum:
         KinematicValue = dunendgarmcSamples[iSample].rw_reco_pi_pMag[iEvent];
         break;
   case kPiTrueMomentum:
         KinematicValue = dunendgarmcSamples[iSample].rw_pi_pMag[iEvent];
         break;
   case kTrueQ2:
	 KinematicValue = dunendgarmcSamples[iSample].rw_Q2[iEvent]; 
	 break;
   case kTrueW:
	 KinematicValue = dunendgarmcSamples[iSample].rw_W[iEvent]; 
	 break;
   case kTrueQ0:
	 KinematicValue = dunendgarmcSamples[iSample].rw_Q0[iEvent]; 
	 break;
   case kTrueQ3:
	 KinematicValue = dunendgarmcSamples[iSample].rw_Q3[iEvent]; 
	 break;
   case kDeltaRecoEnergyThreshold:
         KinematicValue = dunendgarmcSamples[iSample].rw_etrurec_nopionthreshold[iEvent] - dunendgarmcSamples[iSample].rw_etrurec[iEvent];
         break;
   case kIsAccepted:
	 KinematicValue = dunendgarmcSamples[iSample].is_accepted[iEvent]; 
	 break;
   case kIsCC:
         KinematicValue = dunendgarmcSamples[iSample].rw_isCC[iEvent];
         break;
   case kPiTrueMinEnergy:
         KinematicValue = dunendgarmcSamples[iSample].rw_pi_min_energy[iEvent];
         break;
   case kMomResNonAccepted:
	 KinematicValue = dunendgarmcSamples[iSample].momres_nonaccepted[iEvent]; 
	 break;
   case kPDGNonAccepted:
	 KinematicValue = dunendgarmcSamples[iSample].pdg_nonaccepted[iEvent]; 
	 break;
   case kMuonZAngle:
         KinematicValue = dunendgarmcSamples[iSample].muon_z_angle[iEvent];
         break;
   case kRejectedParticleThetaAngle:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_theta_angle[iEvent];
         break;
   case kRejectedParticleMomentum:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_ptot[iEvent];
         break;
   case kRejectedParticleTransverseMomentum:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_pT[iEvent];
         break;
   case kRejectedParticleRadCurvature:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_radcurvature[iEvent];
         break;
   case kSigmaMom:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_sigmamom[iEvent];
         break;
   case kSigmaTheta:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_sigmatheta[iEvent];
         break;
   case kHighestpTParticleTransverseMomentum:
         KinematicValue = dunendgarmcSamples[iSample].highestpart_pT[iEvent];
         break;
   case kHighestpTParticleThetaAngle:
         KinematicValue = dunendgarmcSamples[iSample].highestpart_theta_angle[iEvent];
         break;
   case kHighestpTLengthTrackX:
         KinematicValue = dunendgarmcSamples[iSample].highestpart_lengthtrackx[iEvent];
         break;
   case kHighestpTLengthTrackYZ:
         KinematicValue = dunendgarmcSamples[iSample].highestpart_lengthtrackyz[iEvent];
         break;
   case kRejectedParticleTrackThetaAngle:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_track_theta_angle[iEvent];
         break;
   case kRejectedParticleRatioRadCurvature:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_ratioradcurvature[iEvent];
         break;
   case kTrueSquaredRad:
         KinematicValue = pow(dunendgarmcSamples[iSample].rw_rad[iEvent], 2);
         break;
   case kRejectedParticleBeta:
         KinematicValue = dunendgarmcSamples[iSample].rejectedpart_beta[iEvent];
         break;
   case kECalDepositFrac:
         KinematicValue = dunendgarmcSamples[iSample].ecaldepositfraction[iEvent];
         break;
   case kParticlePDG:
         KinematicValue = dunendgarmcSamples[iSample].particlepdg[iEvent];
         break;
   case kParticleEnergy:
         KinematicValue = dunendgarmcSamples[iSample].particleenergy[iEvent];
         break;
   case kParticleMom:
         KinematicValue = dunendgarmcSamples[iSample].particlemomentum[iEvent];
         break;
   case kParticleDeDx:
         KinematicValue = dunendgarmcSamples[iSample].particlededx[iEvent];
         break;
   case kParticleEvent:
         KinematicValue = dunendgarmcSamples[iSample].particleevent[iEvent];
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
    case kRecoNeutrinoEnergy:
    case kTrueNeutrinoEnergy:
    case kIdealNeutrinoRecoEnergy:
/*         for(double ibins =0; ibins<10*10; ibins++){
           double binval = ibins/10;
           binningVector.push_back(binval);
         }*/
         binningVector = {0.,  0.5,  1.,  1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 5., 6., 10.};
	 break;
    case kTrueQ2:
    case kTrueW:
         for(double ibins =0; ibins<10*50; ibins++){
           double binval = ibins/50;
           binningVector.push_back(binval);
         }
	 break;
    case kTrueQ0:
    case kTrueQ3:
         for(double ibins =0; ibins<3*50+1; ibins++){
           double binval = ibins/50;
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
           binningVector.push_back(ibins-277+TPC_centre_y);
         }
	 break;
    case kRecoZPos:
    case kTrueZPos:
 	 for(double ibins =0; ibins<277*2; ibins++){
           binningVector.push_back(ibins-277+TPC_centre_z);
         }
	 break;
    case kNChargedPions:
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
   case kIsnumuCC: 
   case kInFDV:
         for(double ibins =0; ibins<3; ibins++){
           binningVector.push_back(ibins);
         }
         break;
   case kNMuonsRecoOverTruth:
   case kTrueMinusIdealRecoEnergyRatio:
         for(double ibins =0; ibins<2*100; ibins++){
           binningVector.push_back(-1+(double)(ibins)/100);
         }
         break;
   case kTrueMinusRecoEnergyRatio:
         for(double ibins =0; ibins<20*10; ibins++){
           binningVector.push_back(-10+(double)(ibins)/10);
         }
         break;
   case kTrueMinusIdealRecoEnergy:
   case kDeltaRecoEnergyThreshold:
         for(double ibins =0; ibins<3*100; ibins++){
           binningVector.push_back(-1+(double)(ibins)/100);
         }
         break;
   case kTrueMinusRecoEnergy:
         for(double ibins =0; ibins<20*10; ibins++){
           binningVector.push_back(-10+(double)(ibins)/10);
         }
         break;
   case kNTruePrimMuons:
   case kNTrueMuons:
   case kNRecoMuons:
   case kNRecoPions:
         for(double ibins =0; ibins<10; ibins++){
           binningVector.push_back(ibins);
         }
         break;
    case kRecoLepEnergy:
        for(double ibins =0; ibins<10*10; ibins++){
           binningVector.push_back((double)(ibins)/10);
         } 
         break;
//    case kTrueNeutrinoEnergy:
    case kParticleMom:
        for(double ibins =0; ibins<50*3; ibins++){
           binningVector.push_back((double)pow(10, -2 + ibins/50));
         } 
         break;
    case kParticleEnergy:  
        for(double ibins =0; ibins<20*10; ibins++){
           binningVector.push_back((double)(ibins)/20);
         } 
         break;
    case kTrueLepEnergy:
    case kHighestpTParticleTransverseMomentum:
    case kRejectedParticleTransverseMomentum:
    case kRejectedParticleMomentum:
        for(double ibins =0; ibins<20*10; ibins++){
           binningVector.push_back((double)(ibins)/10);
         } 
         break;
    case kRejectedParticleBeta:
        for(double ibins =0; ibins<1*20+0.05; ibins++){
           binningVector.push_back((double)(ibins)/20);
         } 
         break;
    case kRejectedParticleRadCurvature:
         for(double ibins =0; ibins<50; ibins++){
           binningVector.push_back(ibins);
         }
         break;
    case kRejectedParticleRatioRadCurvature:
         for(double ibins =0; ibins<2*50; ibins++){
           binningVector.push_back(ibins/50);
         }
         break;
    case kHighestpTLengthTrackX:
    case kHighestpTLengthTrackYZ:
        for(double ibins =0; ibins<7*50; ibins++){
           binningVector.push_back(ibins/50);
        }
        break;
    case kTrueRad:
    case kRecoRad:
        for(double ibins =0; ibins<300; ibins = ibins+5){
           binningVector.push_back(ibins);
        }
        break;
    case kTrueSquaredRad:
        for(double ibins =0; ibins<TPCFidRadius*TPCFidRadius + 100; ibins=ibins+100){
           binningVector.push_back(ibins);
        }
        break;
    case kSigmaMom:
        for(double ibins =0; ibins<2*100; ibins++){
           binningVector.push_back((double)(ibins)/100);
        }
        break;
    case kLepRecoPT:
    case kLepRecoPZ:
    case kLepPT:
    case kLepPZ:
        for(double ibins =0; ibins<10*10; ibins++){
           binningVector.push_back((double)(ibins)/10);
        }
        break;
    case kSigmaTheta:
        for(double ibins =0; ibins<2*100; ibins++){
           binningVector.push_back((double)(ibins)/100);
        }
        break;
    case kRejectedParticleTrackThetaAngle:
/*        for(double ibins =0; ibins<20*50; ibins++){
           binningVector.push_back((double)(ibins)/50);
        }
        break;*/
    case kHighestpTParticleThetaAngle:
    case kRejectedParticleThetaAngle:
//    case kRejectedParticleTrackThetaAngle:
/*        for(double ibins = 0; ibins<2*20+1;ibins++){
          binningVector.push_back((double)(ibins/20)-1);
        }
        break;*/
    case kMuonZAngle:
    case kPiZAngle:
    case kPiZRecoAngle:
    case kMuonPiAngle:
    case kMuonPiRecoAngle:
        for(double ibins =0; ibins<360; ibins = ibins+5){
           binningVector.push_back((double)(ibins));
        }
        break;
    case kPiTrueMomentum:
    case kPiRecoMomentum:
    case kPiTrueEnergy:
    case kPiRecoEnergy:
        binningVector = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.75, 1.0};
//        for(double ibins =0; ibins<1.0*100; ibins++){
//           if(ibins>=10){binningVector.push_back((double)(ibins/100)); ibins = ibins+2;}
//           if(ibins>=20){binningVector.push_back((double)(ibins/100)); ibins = ibins+5;}
//           else{binningVector.push_back((double)(ibins/100));}
//        }
        break;
    case kPiTrueMinEnergy:
        binningVector = {0.0, 0.001, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.75, 1.0};
//        for(double ibins =0; ibins<0.1*100; ibins++){
//           binningVector.push_back((double)(ibins/100));
//        }

        break;
    case kMomResNonAccepted:
        binningVector = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
        break;
    case kParticlePDG:
    case kPDGNonAccepted:
        for(double ibins =0; ibins<5000; ibins = ibins+100){
           binningVector.push_back(ibins-2500);
        }
        break;
    case kParticleDeDx:
        for(double ibins =0; ibins<2*100; ibins++){
           binningVector.push_back(ibins/1000);
        }
        break;

    case kECalDepositFrac:
        for(double ibins =0; ibins<1*100; ibins++){
           binningVector.push_back(ibins/100);
        }
        break;
    default:
         for(double ibins =0; ibins<10*5; ibins++){
           binningVector.push_back(ibins/5);
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
//    std::cout<<"StoredSelection: "<<StoredSelection[iStoredSelection][0]<<std::endl;
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
  if(KinematicVar1 != "ECalDepositFrac" && KinematicVar1 != "ParticlePDG" && KinematicVar1 != "ParticleMom" && KinematicVar1 != "ParticleEnergy" && KinematicVar1 != "ParticleDeDx"){
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
//      std::cout<<"First Weight: "<<Weight<<std::endl;
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

//      std::cout<<"Final Weight: "<<Weight<<std::endl;
	  double Var1_Val;
          KinematicTypes Var1 = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicVar1));
	  if (Var1==kPDFBinning) {
		Var1_Val = *(MCSamples[i].x_var[j]);
	  } else {
		Var1_Val = ReturnKinematicParameter(KinematicVar1,i,j);
	  }
	  if (Var1_Val!=__DEFAULT_RETURN_VAL__) {
//                std::cout<<"filling here"<<std::endl;
//                std::cout<<"Var1_Val: "<<Var1_Val<<" Weight: "<<Weight<<std::endl;
		_h1DVar->Fill(Var1_Val,Weight);
	  }
    }
  }
  }
  else{
  for (int i=0;i<getNMCSamples();i++) {
    if (fChannel && (i!=kChannelToFill)) {
      continue;
    }
    for(int j=0;j<getNEventsInSample(i)*7;j++) {
      int eventnum = ReturnKinematicParameter("ParticleEvent",i,j);
      //DB Determine which events pass selection
      if (!IsEventSelected(SelectionStr,i,eventnum)) {
		continue;
      }
      if(MCSamples[i].isNC[eventnum]) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
        MCSamples[i].osc_w[eventnum] = 1.;
      }

      double Weight = getEventWeight(i,eventnum);
//      std::cout<<"First Weight: "<<Weight<<std::endl;
	  if (WeightStyle==1) {
            Weight = 1.0;
            for (int iParam=0; iParam<MCSamples[i].ntotal_weight_pointers[eventnum] ; ++iParam) {
	      Weight *= *(MCSamples[i].total_weight_pointers[eventnum][iParam]);
            }
	  }
	  //ETA - not sure about this
	  //if (MCSamples[i].xsec_w[j] == 0.) continue;
          //NK changing depending on weight style
          if (WeightStyle != 0 && MCSamples[i].xsec_w[eventnum] == 0.) continue;
          if (WeightStyle == 0){Weight = 1.0;}

//      std::cout<<"Final Weight: "<<Weight<<std::endl;
	  double Var1_Val;
          KinematicTypes Var1 = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicVar1));
	  if (Var1==kPDFBinning) {
		Var1_Val = *(MCSamples[i].x_var[eventnum]);
	  } else {
		Var1_Val = ReturnKinematicParameter(KinematicVar1,i,j);
	  }
	  if (Var1_Val!=__DEFAULT_RETURN_VAL__) {
//                std::cout<<"filling here"<<std::endl;
//                std::cout<<"Var1_Val: "<<Var1_Val<<" Weight: "<<Weight<<std::endl;
		_h1DVar->Fill(Var1_Val,Weight);
	  }
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
  std::cout<<"KinematicVars: "<<KinematicVar1<<" "<<KinematicVar2<<std::endl;
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
//  std::cout<<"before h2dvar"<<std::endl;
  TH2D* _h2DVar;
  int kPDFBinning = kNKinematicParams; //XVar
  std::string modename = "all";
  if(SelectionVec.size() != 0){modename = std::to_string((int)(Selection[0][1]));}
  std::string HistName = KinematicVar1 +"_"+KinematicVar2+"_NDGAr_"+modename;
  std::cout<<"HistName: "<<HistName<<std::endl;
  if (ReturnKinematicParameterFromString(KinematicVar1)==kPDFBinning) {
    _h2DVar = (TH2D*)_hPDF2D->Clone(HistName.c_str());
    _h2DVar->Reset();
//    std::cout<<"here 1"<<std::endl;
  } else {
    if (Axis) {
//      std::cout<<"here 2"<<std::endl;
      _h2DVar = new TH2D(HistName.c_str(), KinematicVar1.c_str(),Axis->GetNbins(),Axis->GetXbins()->GetArray(), Axis2->GetNbins(),Axis2->GetXbins()->GetArray());
//    std::cout<<"here 3"<<std::endl;
    } else {
//      std::cout<<"here 4"<<std::endl;
      std::vector<double> xBinEdges = ReturnKinematicParameterBinning(KinematicVar1);
//      std::cout<<"here 5"<<std::endl;
      std::vector<double> yBinEdges = ReturnKinematicParameterBinning(KinematicVar2);
//      std::cout<<"here 6"<<std::endl;
      _h2DVar = new TH2D(HistName.c_str(), KinematicVar1.c_str(), xBinEdges.size()-1, xBinEdges.data(), yBinEdges.size()-1, yBinEdges.data());
//      std::cout<<"here 7"<<std::endl;
    }
  }
//  std::cout<<"made h2DVar"<<std::endl;
  _h2DVar->GetXaxis()->SetTitle(KinematicVar1.c_str());
  _h2DVar->GetYaxis()->SetTitle(KinematicVar2.c_str());
   //This should be the same as FillArray in core basically, except that
  //events will end up in different bins
  for (int i=0;i<getNMCSamples();i++) {
    if (fChannel && (i!=kChannelToFill)) {
      continue;
    }
//    std::cout<<"i "<<i<<std::endl;
    int jmax = getNEventsInSample(i);
//    std::cout<<"jmax "<<jmax<<std::endl;
    bool perparticle = false;
    bool perparticlekin1 = false;
    bool perparticlekin2 = false;
   if(KinematicVar1 == "ECalDepositFrac" || KinematicVar1 == "ParticlePDG" || KinematicVar1 == "ParticleMom" || KinematicVar1 == "ParticleEnergy" || KinematicVar1 == "ParticleDeDx"){jmax = getNEventsInSample(i)*7; perparticle=true; perparticlekin1 = true;}
   if(KinematicVar2 == "ECalDepositFrac" ||  KinematicVar2 == "ParticlePDG" || KinematicVar2 == "ParticleMom" || KinematicVar2 == "ParticleEnergy" || KinematicVar2 == "ParticleDeDx"){jmax = getNEventsInSample(i)*7; perparticle=true; perparticlekin2 = true;}

    for(int j=0;j<jmax;j++) {
      int eventnum = j;
      if(perparticle){eventnum = ReturnKinematicParameter("ParticleEvent",i,j);}
//      std::cout<<"eventnum: "<<eventnum<<std::endl;
      //DB Determine which events pass selection
        if (!IsEventSelected(SelectionStr,i,eventnum)) {
  		continue;
        }
        if(MCSamples[i].isNC[eventnum]) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
          MCSamples[i].osc_w[eventnum] = 1.;
        }

        double Weight = getEventWeight(i,eventnum);
	  if (WeightStyle==1) {
            Weight = 1.0;
            for (int iParam=0; iParam<MCSamples[i].ntotal_weight_pointers[eventnum] ; ++iParam) {
              Weight *= *(MCSamples[i].total_weight_pointers[eventnum][iParam]);
            }
          }
	  //ETA - not sure about this
    	  //if (MCSamples[i].xsec_w[j] == 0.) continue;
          //NK changing depending on weight style
          if (WeightStyle != 0 && MCSamples[i].xsec_w[eventnum] == 0.) continue;
          if (WeightStyle == 0){Weight = 1.0;}
          int eventnum1 = j;
          if(!perparticlekin1 && perparticle){eventnum1 = ReturnKinematicParameter("ParticleEvent",i,j);}
          int eventnum2 = j;
          if(!perparticlekin2 && perparticle){eventnum2 = ReturnKinematicParameter("ParticleEvent",i,j);}
	  double Var1_Val;
          double Var2_Val;
          KinematicTypes Var1 = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicVar1));
          KinematicTypes Var2 = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicVar2));
	  if (Var1==kPDFBinning) {
		Var1_Val = *(MCSamples[i].x_var[eventnum1]);
	  } else {
		Var1_Val = ReturnKinematicParameter(KinematicVar1,i,eventnum1);
	  }
	  if (Var2==kPDFBinning) {
		Var2_Val = *(MCSamples[i].y_var[eventnum2]);
	  } else {
		Var2_Val = ReturnKinematicParameter(KinematicVar2,i,eventnum2);
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
