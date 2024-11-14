#include "samplePDFDUNEBeamFD.h"

#include "manager/manager.h"

#include "splines/splinesDUNE.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

samplePDFDUNEBeamFD::samplePDFDUNEBeamFD(std::string mc_version_,
                                         covarianceXsec *xsec_cov_)
    : samplePDFFDBase(mc_version_, xsec_cov_) {
  // Call insitialise in samplePDFFD
  Initialise();
}

samplePDFDUNEBeamFD::~samplePDFDUNEBeamFD() {}

void samplePDFDUNEBeamFD::Init() {
  dunemcSamples.resize(nSamples, dunemc_base());

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");
}

void samplePDFDUNEBeamFD::SetupSplines() {

  ///@todo move all of the spline setup into core
  if (XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0) {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will create a spline object",
        XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splinesDUNE *DUNESplines = new splinesDUNE(XsecCov);
    splineFile = (splineFDBase *)DUNESplines;
    InitialiseSplineObject();
  } else {
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or "
                  "evaluate splines",
                  XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splineFile = nullptr;
  }
}

void samplePDFDUNEBeamFD::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    for (int ev_i = 0; ev_i < dunemcSamples[i].nEvents; ++ev_i) {
      MCSamples[i].ntotal_weight_pointers[ev_i] = 6;
      MCSamples[i].total_weight_pointers[ev_i] =
          new const double *[MCSamples[i].ntotal_weight_pointers[ev_i]];
      MCSamples[i].total_weight_pointers[ev_i][0] = &(dunemcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[ev_i][1] = &(dunemcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[ev_i][2] =
          MCSamples[i].osc_w_pointer[ev_i];
      MCSamples[i].total_weight_pointers[ev_i][3] =
          &(dunemcSamples[i].rw_berpaacvwgt[ev_i]);
      MCSamples[i].total_weight_pointers[ev_i][4] =
          &(dunemcSamples[i].flux_w[ev_i]);
      MCSamples[i].total_weight_pointers[ev_i][5] =
          &(MCSamples[i].xsec_w[ev_i]);
    }
  }
}

int samplePDFDUNEBeamFD::setupExperimentMC(int iSample) {

  dunemc_base *duneobj = &(dunemcSamples[iSample]);

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files[iSample]);

  std::unique_ptr<TFile> in_file{
      TFile::Open(mc_files[iSample].c_str(), "READ")};
  if (!in_file) {
    MACH3LOG_ERROR("Could not open file: {}", mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  auto in_tree = in_file->Get<TTree>("cafmaker/cafTree");
  if (!in_tree) {
    MACH3LOG_ERROR("Could not find \"cafmaker/cafTree\" TTree in {}",
                   mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TTreeReader tree_rdr(in_tree);
  // The branch "px" contains floats; access them as myPx.
  TTreeReaderValue<caf::StandardRecord> sr(tree_rdr, "rec");

  duneobj->nEvents = tree_rdr.GetEntries();

  MACH3LOG_INFO("Found \"cafmaker/cafTree\" TTree in {}", mc_files[iSample]);
  MACH3LOG_INFO("With number of entries: {}", duneobj->nEvents);

  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu_shifted = new double[duneobj->nEvents];
  duneobj->rw_cvnnue_shifted = new double[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];

  duneobj->nutype = sample_nutype[iSample];
  duneobj->oscnutype = sample_oscnutype[iSample];
  duneobj->signal = sample_signal[iSample];

  // after allocation, set all event entries to a bad bad value
  std::fill_n(duneobj->rw_etru, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_erec_shifted, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_vtx_x, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_vtx_y, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_vtx_z, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_cvnnumu_shifted, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_cvnnue_shifted, duneobj->nEvents, -99999);
  std::fill_n(duneobj->mode, duneobj->nEvents, -99999);
  std::fill_n(duneobj->Target, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_isCC, duneobj->nEvents, -99999);

  for (auto entryi : tree_rdr) {

    if (sr->mc.nu.size() < 1) { // no true neutrinos here
      continue;
    }

    duneobj->rw_etru[entryi] = sr->mc.nu[0].E;
    duneobj->mode[entryi] = sr->mc.nu[0].mode;
    duneobj->Target[entryi] = (sr->mc.nu[0].targetPDG / 1000) % 1000;
    duneobj->rw_isCC[entryi] = sr->mc.nu[0].iscc;

    if (sr->common.ixn.pandora.size() < 1) { // no reconstructed objects
      continue;
    }

    duneobj->rw_erec_shifted[entryi] = sr->common.ixn.pandora[0].Enu.lep_calo;
    duneobj->rw_vtx_x[entryi] = sr->common.ixn.pandora[0].vtx.x;
    duneobj->rw_vtx_y[entryi] = sr->common.ixn.pandora[0].vtx.y;
    duneobj->rw_vtx_z[entryi] = sr->common.ixn.pandora[0].vtx.z;
    duneobj->rw_cvnnumu_shifted[entryi] =
        sr->common.ixn.pandora[0].nuhyp.cvn.numu;
    duneobj->rw_cvnnue_shifted[entryi] =
        sr->common.ixn.pandora[0].nuhyp.cvn.nue;
  }

  return duneobj->nEvents;
}

TH1D *samplePDFDUNEBeamFD::get1DVarHist(KinematicTypes Var1, int kModeToFill,
                                        int kChannelToFill, int WeightStyle,
                                        TAxis *Axis) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill != -1) {
    if (kChannelToFill > dunemcSamples.size()) {
      MACH3LOG_ERROR("Required channel is not available. kChannelToFill should "
                     "be between 0 and {}",
                     dunemcSamples.size());
      MACH3LOG_ERROR("kChannelToFill given: {}", kChannelToFill);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill != -1) {
    if (kModeToFill > kMaCh3_nModes) {
      MACH3LOG_ERROR("Required mode is not available. kModeToFill should be "
                     "between 0 and {}",
                     kMaCh3_nModes);
      MACH3LOG_ERROR("kModeToFill given: {}", kModeToFill);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fMode = true;
  } else {
    fMode = false;
  }

  std::vector<std::vector<double>> SelectionVec;

  if (fMode) {
    std::vector<double> SelecMode(3);
    SelecMode[0] = kM3Mode;
    SelecMode[1] = kModeToFill;
    SelecMode[2] = kModeToFill + 1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    std::vector<double> SelecChannel(3);
    SelecChannel[0] = kOscChannel;
    SelecChannel[1] = kChannelToFill;
    SelecChannel[2] = kChannelToFill + 1;
    SelectionVec.push_back(SelecChannel);
  }

  return get1DVarHist(Var1, SelectionVec, WeightStyle, Axis);
}

/*! DB New version of get1DVarHist which only fills histogram with events
 * passing IsEventSelected This works by having the Selection vector, where each
 * component of Selection is a 2 or 3 length vector If Selection[i].size()==3,
 * Selection[i][0] is the ND280KinematicType which is being cut, and only events
 * with ND280KinematicType values between Selection[i][1] and Selection[i][2]
 * are accepted
 */
TH1D *
samplePDFDUNEBeamFD::get1DVarHist(KinematicTypes Var1,
                                  std::vector<std::vector<double>> SelectionVec,
                                  int WeightStyle, TAxis *Axis) {

  Selection = SelectionVec;

  for (unsigned int iStoredSelection = 0;
       iStoredSelection < StoredSelection.size(); iStoredSelection++) {
    Selection.push_back(StoredSelection[iStoredSelection]);
  }

  for (unsigned int iSelection = 0; iSelection < Selection.size();
       iSelection++) {
    if (Selection[iSelection].size() != 3) {
      MACH3LOG_ERROR("Selection Vector[{}] is not formed correctly. Expect "
                     "size == 3, given: {}",
                     iSelection, Selection[iSelection].size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  // DB Cut on OscChannel in this function due to speed increase from
  // considering duneSamples structure (ie. Array of length NChannels)
  bool fChannel = false;
  int kChannelToFill = -1;
  for (unsigned int iSelection = 0; iSelection < Selection.size();
       iSelection++) {
    if (Selection[iSelection][0] == kOscChannel) {
      fChannel = true;
      kChannelToFill = Selection[iSelection][1];
    }
  }

  if (fChannel && kChannelToFill > dunemcSamples.size()) {
    MACH3LOG_ERROR("Required channel is not available. kChannelToFill should "
                   "be between 0 and {}",
                   dunemcSamples.size());
    MACH3LOG_ERROR("kChannelToFill given: {}", kChannelToFill);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TH1D *_h1DVar;
  std::vector<double> xBinEdges =
      ReturnKinematicParameterBinning(ReturnStringFromKinematicParameter(Var1));
  _h1DVar = new TH1D("", "", xBinEdges.size() - 1, xBinEdges.data());

  // This should be the same as FillArray in core basically, except that
  // events will end up in different bins
  for (int i = 0; i < dunemcSamples.size(); i++) {
    if (fChannel && (i != kChannelToFill)) {
      continue;
    }
    for (int ev_i = 0; ev_i < dunemcSamples[i].nEvents; ev_i++) {

      // DB Determine which events pass selection
      if (!IsEventSelected(i, ev_i)) {
        continue;
      }

      double Weight = GetEventWeight(i, ev_i);
      if (WeightStyle == 1) {
        Weight = *(MCSamples[i].osc_w_pointer[ev_i]) * dunemcSamples[i].pot_s *
                 dunemcSamples[i].norm_s * dunemcSamples[i].flux_w[ev_i];
      }

      // ETA - not sure about this
      if (MCSamples[i].xsec_w[ev_i] == 0.)
        continue;

      double Var1_Val;

      Var1_Val = ReturnKinematicParameter(Var1, i, ev_i);
      if (Var1_Val != _DEFAULT_RETURN_VAL_) {
        _h1DVar->Fill(Var1_Val, Weight);
      }
    }
  }

  /* DB: This is commented out be default
  // This code shifts the histogram meaning to Events/Bin Width but this affects
  the overall integral of the histogram so it should not be used anywhere we
  care about event rates
  // We could use Hist->Integral("width") but it would require a lot of
  modification throughout the code

  if (Var1!=kPDFBinning) {
    //_h1DVar->SetBinContent(1,_h1DVar->GetBinContent(0)+_h1DVar->GetBinContent(1));
    //_h1DVar->SetBinContent(_h1DVar->GetNbinsX(),_h1DVar->GetBinContent(_h1DVar->GetNbinsX())+_h1DVar->GetBinContent(_h1DVar->GetNbinsX()+1));

    for (int x=1;x<=_h1DVar->GetNbinsX();x++) {
      _h1DVar->SetBinContent(x,_h1DVar->GetBinContent(x)/_h1DVar->GetXaxis()->GetBinWidth(x));
    }

    _h1DVar->GetYaxis()->SetTitle("Events/Bin Width");
  }
  */

  return _h1DVar;
}

double samplePDFDUNEBeamFD::ReturnKinematicParameter(double KinematicVariable,
                                                     int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicVariable, iSample, iEvent);
}

double
samplePDFDUNEBeamFD::ReturnKinematicParameter(std::string KinematicParameter,
                                              int iSample, int iEvent) {
  return ReturnKinematicParameter(
      ReturnKinematicParameterFromString(KinematicParameter), iSample, iEvent);
}

const double *samplePDFDUNEBeamFD::GetPointerToKinematicParameter(
    std::string KinematicParameter, int iSample, int iEvent) {
  return GetPointerToKinematicParameter(
      ReturnKinematicParameterFromString(KinematicParameter), iSample, iEvent);
}

int samplePDFDUNEBeamFD::ReturnKinematicParameterFromString(
    std::string KinematicParameterStr) {
  if (KinematicParameterStr.find("TrueNeutrinoEnergy") != std::string::npos) {
    return kTrueNeutrinoEnergy;
  }
  if (KinematicParameterStr.find("RecoNeutrinoEnergy") != std::string::npos) {
    return kRecoNeutrinoEnergy;
  }
  if (KinematicParameterStr.find("TrueXPos") != std::string::npos) {
    return kTrueXPos;
  }
  if (KinematicParameterStr.find("TrueYPos") != std::string::npos) {
    return kTrueYPos;
  }
  if (KinematicParameterStr.find("TrueZPos") != std::string::npos) {
    return kTrueZPos;
  }
  if (KinematicParameterStr.find("CVNNumu") != std::string::npos) {
    return kCVNNumu;
  }
  if (KinematicParameterStr.find("CVNNue") != std::string::npos) {
    return kCVNNue;
  }
  if (KinematicParameterStr.find("M3Mode") != std::string::npos) {
    return kM3Mode;
  }
}

const double *
samplePDFDUNEBeamFD::GetPointerToKinematicParameter(double KinematicVariable,
                                                    int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes)std::round(KinematicVariable);

  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    return &dunemcSamples[iSample].rw_etru[iEvent];
  case kRecoNeutrinoEnergy:
    return &dunemcSamples[iSample].rw_erec_shifted[iEvent];
  case kTrueXPos:
    return &dunemcSamples[iSample].rw_vtx_x[iEvent];
  case kTrueYPos:
    return &dunemcSamples[iSample].rw_vtx_y[iEvent];
  case kTrueZPos:
    return &dunemcSamples[iSample].rw_vtx_z[iEvent];
  case kCVNNumu:
    return &dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
  case kCVNNue:
    return &dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

inline std::string samplePDFDUNEBeamFD::ReturnStringFromKinematicParameter(
    int KinematicParameter) {
  std::string KinematicString = "";

  switch (KinematicParameter) {
  case kRecoNeutrinoEnergy:
    KinematicString = "RecoNeutrinoEnergy";
    break;
  case kTrueNeutrinoEnergy:
    KinematicString = "RecoNeutrinoEnergy";
    break;
  case kTrueXPos:
    KinematicString = "TrueXPos";
    break;
  case kTrueYPos:
    KinematicString = "TrueYPos";
    break;
  case kTrueZPos:
    KinematicString = "TrueZPos";
    break;
  case kCVNNumu:
    KinematicString = "CVNNumu";
    break;
  case kCVNNue:
    KinematicString = "CVNNue";
    break;
  case kM3Mode:
    KinematicString = "M3Mode";
    break;
  default:
    break;
  }

  return KinematicString;
}

void samplePDFDUNEBeamFD::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  fdmc_base *fdobj = &(MCSamples[iSample]);

  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->SampleDetID = SampleDetID;

  for (int iEvent = 0; iEvent < fdobj->nEvents; ++iEvent) {
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]);
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
  }
}

void samplePDFDUNEBeamFD::applyShifts(int iSample, int iEvent) {}

std::vector<double> samplePDFDUNEBeamFD::ReturnKinematicParameterBinning(
    std::string KinematicParameterStr) {
  std::vector<double> binningVector;
  KinematicTypes KinematicParameter = static_cast<KinematicTypes>(
      ReturnKinematicParameterFromString(KinematicParameterStr));

  int nBins = 0;
  double bin_width = 0;
  switch (KinematicParameter) {
  case (kRecoNeutrinoEnergy):
    nBins = 20;
    bin_width = 0.5; // GeV
    break;
  case (kTrueNeutrinoEnergy):
    nBins = 20;
    bin_width = 0.5; // GeV
    break;
  default:
    nBins = 10;
    bin_width = 1.0;
    break;
  }

  for (int bin_i = 0; bin_i < nBins; bin_i++) {
    binningVector.push_back(bin_i * bin_width);
  }

  return binningVector;
}
