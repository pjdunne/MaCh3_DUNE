#include "samplePDFDUNEBeamFD.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

samplePDFDUNEBeamFD::samplePDFDUNEBeamFD(std::string mc_version_,
                                         covarianceXsec *xsec_cov_)
    : samplePDFFDBase(mc_version_, xsec_cov_) {
  // Call insitialise in samplePDFFD
  Initialise();
}

samplePDFDUNEBeamFD::~samplePDFDUNEBeamFD() {}

void samplePDFDUNEBeamFD::Init() {

  auto const &configyml = SampleManager->raw();

  if (!configyml["DUNESampleBools"] ||
      !configyml["DUNESampleBools"]["iselike"]) {
    MACH3LOG_ERROR(
        "Did not find DUNESampleBools:iselike in {}, please add this",
        SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  iselike = configyml["DUNESampleBools"]["iselike"].as<bool>();

  if (!configyml["POT"]) {
    MACH3LOG_ERROR("POT not defined in {}, please add this!",
                   SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  pot = configyml["POT"].as<double>();

  // create dunemc storage
  dunemcSamples.resize(nSamples);

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");
}

void samplePDFDUNEBeamFD::SetupSplines() {

  ///@todo move all of the spline setup into core
  if (spline_files.size() &&
      (XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0)) {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will create a spline object",
        XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splineFile = new splinesDUNE(XsecCov);
    InitialiseSplineObject();
  } else {
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or "
                  "evaluate splines",
                  XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splineFile = nullptr;
  }

  return;
}

void samplePDFDUNEBeamFD::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      MCSamples[i].total_weight_pointers[j] =
          new const double *[MCSamples[i].ntotal_weight_pointers[j]];
      MCSamples[i].total_weight_pointers[j][0] = &(dunemcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunemcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] =
          &(dunemcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(dunemcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

int samplePDFDUNEBeamFD::setupExperimentMC(int iSample) {

  auto &duneobj = dunemcSamples[iSample];

  int nutype = sample_nutype[iSample];
  int oscnutype = sample_oscnutype[iSample];
  bool signal = sample_signal[iSample];

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files[iSample].native());

  std::unique_ptr<TFile> infile(TFile::Open(mc_files[iSample].c_str()));
  auto tree = infile->Get<TTree>("cafmaker/cafTree");
  if (!tree) {
    MACH3LOG_ERROR("Failed to read cafmaker/cafTree from {} ",
                   mc_files[iSample].native());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  caf::Proxy<caf::StandardRecord> srp(tree,"rec");

  // now fill the actual variables
  duneobj.norm_s = 1;
  duneobj.pot_s = 1; // read pot out of meta tree

  // duneobj.nEvents = _data->GetEntries();
  duneobj.nutype = nutype;
  duneobj.oscnutype = oscnutype;
  duneobj.signal = signal;

  // allocate memory for dunemc variables
  duneobj.rw_cvnnumu.resize(duneobj.nEvents);
  duneobj.rw_cvnnue.resize(duneobj.nEvents);
  duneobj.rw_cvnnumu_shifted.resize(duneobj.nEvents);
  duneobj.rw_cvnnue_shifted.resize(duneobj.nEvents);
  duneobj.rw_etru.resize(duneobj.nEvents);
  duneobj.rw_erec.resize(duneobj.nEvents);
  duneobj.rw_erec_shifted.resize(duneobj.nEvents);
  duneobj.rw_erec_had.resize(duneobj.nEvents);
  duneobj.rw_erec_lep.resize(duneobj.nEvents);

  duneobj.true_q0.resize(duneobj.nEvents);
  duneobj.true_q3.resize(duneobj.nEvents);

  duneobj.rw_eRecoP.resize(duneobj.nEvents);
  duneobj.rw_eRecoPip.resize(duneobj.nEvents);
  duneobj.rw_eRecoPim.resize(duneobj.nEvents);
  duneobj.rw_eRecoPi0.resize(duneobj.nEvents);
  duneobj.rw_eRecoN.resize(duneobj.nEvents);

  duneobj.rw_LepE.resize(duneobj.nEvents);
  duneobj.rw_eP.resize(duneobj.nEvents);
  duneobj.rw_ePip.resize(duneobj.nEvents);
  duneobj.rw_ePim.resize(duneobj.nEvents);
  duneobj.rw_ePi0.resize(duneobj.nEvents);
  duneobj.rw_eN.resize(duneobj.nEvents);

  duneobj.rw_theta.resize(duneobj.nEvents);
  duneobj.flux_w.resize(duneobj.nEvents);
  duneobj.rw_isCC.resize(duneobj.nEvents);
  duneobj.rw_nuPDGunosc.resize(duneobj.nEvents);
  duneobj.rw_nuPDG.resize(duneobj.nEvents);
  duneobj.rw_berpaacvwgt.resize(duneobj.nEvents);
  duneobj.rw_vtx_x.resize(duneobj.nEvents);
  duneobj.rw_vtx_y.resize(duneobj.nEvents);
  duneobj.rw_vtx_z.resize(duneobj.nEvents);

  duneobj.global_bin_number.resize(duneobj.nEvents);

  duneobj.mode.resize(duneobj.nEvents);
  duneobj.Target.resize(duneobj.nEvents);

  bool need_global_bin_numbers = (XVarStr == "global_bin_number");

  // FILL DUNE STRUCT
  for (int i = 0; i < duneobj.nEvents; ++i) { // Loop through tree

    // duneobj.rw_cvnnumu[i] = _cvnnumu;
    // duneobj.rw_cvnnue[i] = _cvnnue;
    // duneobj.rw_cvnnumu_shifted[i] = _cvnnumu;
    // duneobj.rw_cvnnue_shifted[i] = _cvnnue;

    // if (iselike) {
    //   duneobj.rw_erec[i] = _erec_nue;
    //   duneobj.rw_erec_shifted[i] = _erec_nue;
    //   duneobj.rw_erec_had[i] = _erec_had_nue;
    //   duneobj.rw_erec_lep[i] = _erec_lep_nue;
    // } else {
    //   duneobj.rw_erec[i] = _erec;
    //   duneobj.rw_erec_shifted[i] = _erec;
    //   duneobj.rw_erec_had[i] = _erec_had;
    //   duneobj.rw_erec_lep[i] = _erec_lep;
    // }

    // duneobj.true_q0[i] = _ev - _LepE;
    // duneobj.true_q3[i] = (TVector3{_NuMomX, _NuMomY, _NuMomZ} -
    //                       TVector3{_LepMomX, _LepMomY, _LepMomZ})
    //                          .Mag();

    // duneobj.rw_eRecoP[i] = _eRecoP;
    // duneobj.rw_eRecoPip[i] = _eRecoPip;
    // duneobj.rw_eRecoPim[i] = _eRecoPim;
    // duneobj.rw_eRecoPi0[i] = _eRecoPi0;
    // duneobj.rw_eRecoN[i] = _eRecoN;

    // duneobj.rw_LepE[i] = _LepE;
    // duneobj.rw_eP[i] = _eP;
    // duneobj.rw_ePip[i] = _ePip;
    // duneobj.rw_ePim[i] = _ePim;
    // duneobj.rw_ePi0[i] = _ePi0;
    // duneobj.rw_eN[i] = _eN;

    // duneobj.rw_etru[i] = _ev;
    // duneobj.rw_theta[i] = _LepNuAngle;
    // duneobj.rw_isCC[i] = _isCC;
    // duneobj.rw_nuPDGunosc[i] = _nuPDGunosc;
    // duneobj.rw_nuPDG[i] = _nuPDG;
    // duneobj.rw_berpaacvwgt[i] = _BeRPA_cvwgt;
    // duneobj.rw_vtx_x[i] = _vtx_x;
    // duneobj.rw_vtx_y[i] = _vtx_y;
    // duneobj.rw_vtx_z[i] = _vtx_z;

    if (need_global_bin_numbers) {
      duneobj.global_bin_number[i] =
          GetGenericBinningGlobalBinNumber(iSample, i);
    }
    // Assume everything is on Argon for now....
    duneobj.Target[i] = 40;

    // int mode = TMath::Abs(_mode);
    // duneobj.mode[i] = SIMBMode_ToMaCh3Mode(mode, _isCC);

    duneobj.flux_w[i] = 1.0;
  }

  return duneobj.nEvents;
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
  std::vector<double> xBinEdges = ReturnKinematicParameterBinning(Var1);
  _h1DVar = new TH1D("", "", xBinEdges.size() - 1, xBinEdges.data());

  // This should be the same as FillArray in core basically, except that
  // events will end up in different bins
  for (int i = 0; i < dunemcSamples.size(); i++) {
    if (fChannel && (i != kChannelToFill)) {
      continue;
    }
    for (int j = 0; j < dunemcSamples[i].nEvents; j++) {

      // DB Determine which events pass selection
      if (!IsEventSelected(i, j)) {
        continue;
      }

      double Weight = GetEventWeight(i, j);
      if (WeightStyle == 1) {
        Weight = *(MCSamples[i].osc_w_pointer[j]) * dunemcSamples[i].pot_s *
                 dunemcSamples[i].norm_s * dunemcSamples[i].flux_w[j];
      }

      // ETA - not sure about this
      if (MCSamples[i].xsec_w[j] == 0.)
        continue;

      double Var1_Val;

      Var1_Val = ReturnKinematicParameter(Var1, i, j);
      if (Var1_Val != _DEFAULT_RETURN_VAL_) {
        _h1DVar->Fill(Var1_Val, Weight);
      }
    }
  }

  return _h1DVar;
}

double const &samplePDFDUNEBeamFD::ReturnKinematicParameterByReference(
    int KinematicParameter, int iSample, int iEvent) {

  switch (KinematicParameter) {
  case kTrueNeutrinoEnergy:
    return dunemcSamples[iSample].rw_etru[iEvent];
  case kRecoNeutrinoEnergy:
    return dunemcSamples[iSample].rw_erec_shifted[iEvent];
  case kTrueXPos:
    return dunemcSamples[iSample].rw_vtx_x[iEvent];
  case kTrueYPos:
    return dunemcSamples[iSample].rw_vtx_y[iEvent];
  case kTrueZPos:
    return dunemcSamples[iSample].rw_vtx_z[iEvent];
  case kCVNNumu:
    return dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
  case kCVNNue:
    return dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
  case kGlobalBinNumber:
    return dunemcSamples[iSample].global_bin_number[iEvent];
  case kELepRec: {
    return dunemcSamples[iSample].rw_erec_lep[iEvent];
  }
  case kq0:
    return dunemcSamples[iSample].true_q0[iEvent];
  case kq3:
    return dunemcSamples[iSample].true_q3[iEvent];
  default:

    if (KinematicParameter >= kNDefaultProjections) {
      MACH3LOG_ERROR("ReturnKinematicParameterByReference passed "
                     "KinematicParameter: {}, which appears to refer to user "
                     "projection: {}, but user projections cannot be evaluated "
                     "by reference.",
                     KinematicParameter,
                     KinematicParameter - kNDefaultProjections);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    MACH3LOG_ERROR(
        "ReturnKinematicParameterByReference Did not recognise "
        "Kinematic Parameter type: {}. Is it possibly only available via "
        "ReturnKinematicParameter "
        "(not by reference?) if you need it here, give it storage in "
        "dunemc_base and move it to here.",
        KinematicParameter);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

double samplePDFDUNEBeamFD::ReturnKinematicParameter(int KinematicParameter,
                                                     int iSample, int iEvent) {

  if (KinematicParameter >= kNDefaultProjections) {

    if ((KinematicParameter - kNDefaultProjections) > user_projections.size()) {
      MACH3LOG_ERROR(
          "Passed KinematicParameter: {}, which appears to refer to user "
          "projection: {}, but we only have {} user projections defined.",
          KinematicParameter, KinematicParameter - kNDefaultProjections,
          user_projections.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    return user_projections[KinematicParameter - kNDefaultProjections].proj(
        dunemcSamples[iSample], iEvent);
  }

  switch (KinematicParameter) {
  case kERecQE: {
    constexpr double V = 0;        // 0 binding energy for now
    constexpr double mn = 939.565; // neutron mass
    constexpr double mp = 938.272; // proton mass

    double mN_eff = mn - V;
    double mN_oth = mp;

    if (dunemcSamples[iSample].rw_nuPDGunosc[iEvent] <
        0) { // if anti-neutrino, swap target/out masses
      mN_eff = mp - V;
      mN_oth = mn;
    }

    double el = dunemcSamples[iSample].rw_erec_lep[iEvent];

    // this is funky, but don't be scared, it defines an annonymous function
    // in place that grabs the lepton mass in MeV when given the neutrino PDG
    // and whether the interaction was CC or NC and then immediately calls it.
    // It's basically a generalisation of the ternary operator.
    double ml =
        [](int nupdg, bool isCC) {
          switch (std::abs(nupdg)) {
          case 12: {
            return isCC ? 0.511 : 0;
          }
          case 14: {
            return isCC ? 105.66 : 0;
          }
          case 16: {
            return isCC ? 1777.0 : 0;
          }
          }
        }(dunemcSamples[iSample].rw_nuPDGunosc[iEvent],
          dunemcSamples[iSample].rw_isCC[iEvent]);

    double pl = std::sqrt(el * el - ml * ml); // momentum of lepton

    double rEnu =
        (2 * mN_eff * el - ml * ml + mN_oth * mN_oth - mN_eff * mN_eff) /
        (2 * (mN_eff - el +
              pl * std::cos(dunemcSamples[iSample].rw_theta[iEvent])));

    return rEnu;
  }
  case kEHadRec: {

    return dunemcSamples[iSample].rw_eRecoP[iEvent] +
           dunemcSamples[iSample].rw_eRecoPip[iEvent] +
           dunemcSamples[iSample].rw_eRecoPim[iEvent] +
           dunemcSamples[iSample].rw_eRecoPi0[iEvent] +
           dunemcSamples[iSample].rw_eRecoN[iEvent];
  }
  default: {
    return ReturnKinematicParameterByReference(KinematicParameter, iSample,
                                               iEvent);
  }
  }
}

int samplePDFDUNEBeamFD::ReturnKinematicParameterFromString(
    std::string KinematicParameterStr) {
  if (KinematicParameterStr == "TrueNeutrinoEnergy") {
    return kTrueNeutrinoEnergy;
  }
  if (KinematicParameterStr == "RecoNeutrinoEnergy") {
    return kRecoNeutrinoEnergy;
  }
  if (KinematicParameterStr == "TrueXPos") {
    return kTrueXPos;
  }
  if (KinematicParameterStr == "TrueYPos") {
    return kTrueYPos;
  }
  if (KinematicParameterStr == "TrueZPos") {
    return kTrueZPos;
  }
  if (KinematicParameterStr == "CVNNumu") {
    return kCVNNumu;
  }
  if (KinematicParameterStr == "CVNNue") {
    return kCVNNue;
  }
  if (KinematicParameterStr == "M3Mode") {
    return kM3Mode;
  }
  if (KinematicParameterStr == "global_bin_number") {
    return kGlobalBinNumber;
  }
  if (KinematicParameterStr == "q0") {
    return kq0;
  }
  if (KinematicParameterStr == "q3") {
    return kq3;
  }
  if (KinematicParameterStr == "ERecQE") {
    return kERecQE;
  }
  if (KinematicParameterStr == "ELepRec") {
    return kELepRec;
  }
  if (KinematicParameterStr == "EHadRec") {
    return kEHadRec;
  }

  for (size_t up_it = 0; up_it < user_projections.size(); ++up_it) {
    if (KinematicParameterStr == user_projections[up_it].name) {
      return kNDefaultProjections + up_it;
    }
  }

  std::stringstream ss;
  ss << "[ERROR]: " << __FILE__ << ":" << __LINE__
     << "failed to translate kinematic parameter string "
     << KinematicParameterStr << " to parameter id.";
  throw std::runtime_error(ss.str());
}

std::string samplePDFDUNEBeamFD::ReturnStringFromKinematicParameter(
    int KinematicParameter) {

  switch (KinematicParameter) {
  case kRecoNeutrinoEnergy:
    return "RecoNeutrinoEnergy";
  case kTrueNeutrinoEnergy:
    return "RecoNeutrinoEnergy";
  case kTrueXPos:
    return "TrueXPos";
  case kTrueYPos:
    return "TrueYPos";
  case kTrueZPos:
    return "TrueZPos";
  case kCVNNumu:
    return "CVNNumu";
  case kCVNNue:
    return "CVNNue";
  case kM3Mode:
    return "M3Mode";
  case kGlobalBinNumber:
    return "global_bin_number";
  case kq0:
    return "q0";
  case kq3:
    return "q3";
  case kERecQE:
    return "ERecQE";
  case kELepRec:
    return "ELepRec";
  case kEHadRec:
    return "EHadRec";
  default: {
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  }
}

void samplePDFDUNEBeamFD::setupFDMC(int iSample) {
  auto &duneobj = dunemcSamples[iSample];
  fdmc_base *fdobj = &(MCSamples[iSample]);

  fdobj->nutype = duneobj.nutype;
  fdobj->oscnutype = duneobj.oscnutype;
  fdobj->signal = duneobj.signal;
  fdobj->SampleDetID = SampleDetID;

  for (int iEvent = 0; iEvent < fdobj->nEvents; ++iEvent) {
    fdobj->rw_etru[iEvent] = &(duneobj.rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj.mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj.Target[iEvent]);
    fdobj->isNC[iEvent] = !(duneobj.rw_isCC[iEvent]);
  }
}

void samplePDFDUNEBeamFD::applyShifts(int iSample, int iEvent) {}

std::vector<double>
samplePDFDUNEBeamFD::ReturnKinematicParameterBinning(int KinematicParameter) {
  std::vector<double> binningVector;

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

int samplePDFDUNEBeamFD::AddProjection(
    std::string const &name,
    std::function<double(dunemc_base const &, int)> proj) {
  for (auto const &up : user_projections) {
    if (name == up.name) {
      MACH3LOG_ERROR("Attempting to overwrite existing UserProjection: {}, "
                     "please pick a new name",
                     name);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  user_projections.push_back(UserProjection{name, proj});
  return kNDefaultProjections + (user_projections.size() - 1);
}

std::vector<samplePDFDUNEBeamFD::UserProjection>
    samplePDFDUNEBeamFD::user_projections;
