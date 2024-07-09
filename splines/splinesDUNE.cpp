#include "splinesDUNE.h"
#include "samplePDFDUNE/StructsDUNE.h"

splinesDUNE::splinesDUNE(const char *spline, int nutype, int nevents, int DetID, covarianceXsec* xsec_cov) : splineFDBase(spline, nutype, nevents, DetID, xsec_cov) // constructor for etrue-var1 splines
{

  FindUniqueModes();

  SetupSplines();

}

splinesDUNE::splinesDUNE(const char *spline, int nutype, int nevents, double BinningOpt, int DetID, covarianceXsec* xsec_cov) : splineFDBase(spline, nutype, nevents, BinningOpt, DetID, xsec_cov){
  // constructor for etrue-var1-var2 splines 
  FindUniqueModes();

  SetupSplines(BinningOpt);
}


// ---- SetupSplines (first: original erec version, then 2d version) ---- //
void splinesDUNE::SetupSplines()
{
  // ETA - need to think about how to do this configurably. If we store all the dev_blah_sp in one vector then can loop through giving the address of each?
  // Also need to pass in the name of the spline in the splinefile from the xsec xml to the covarianceXsec class
  std::vector<syst*> systs;

#if USE_SPLINE_FD == USE_TSpline3_red_FD
  std::cout << "###########################" << std::endl;
  std::cout << "USING TSPLINE3 RED !!!!!" << std::endl;
  std::cout << "###########################" << std::endl;
#endif

  // Set spline binning
  SetSplineBinning();

  //vector to keep track of which splines are flat. We use this later on to make
  //sure we've loaded everything correctly 
  std::vector<std::vector<std::vector<std::vector<bool> > > > flat_vec;

  //ETA - testing new Setup
  //DB Now use detid to determine number of spline systematics, names and corresponding modes
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::cout << "Length of SplineFileParsNames is " << SplineFileParsNames.size() << std::endl;
  std::cout << "SampleDetID is " << SampleDetID << std::endl;
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);

  std::cout << "Expecting " << numSplineParams << " splines " << std::endl;

  for(int isyst=0; isyst<numSplineParams ; isyst++){  // loop over systematics 
	//std::cout << "On isyst " << isyst << std::endl;
#if USE_SPLINE_FD == USE_TSpline3_FD
	std::vector<std::vector<std::vector<TSpline3*> > > tmp_tmp_imode;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	std::vector<std::vector<std::vector<TSpline3_red*> > > tmp_tmp_imode;
#endif
	std::vector<std::vector<std::vector<bool> > > tmp_flat_mode; 
	//ETA adding in this to store weights for all splines
	std::vector<std::vector<std::vector<double> > > tmp_w_mode; 
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
#if USE_SPLINE_FD == USE_TSpline3_FD
	  std::vector<std::vector<TSpline3*> > tmp_mbin;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	  std::vector<std::vector<TSpline3_red*> > tmp_mbin;
#endif
	  std::vector<std::vector<bool> > tmp_flat_enu;
	  std::vector<std::vector<double> > tmp_w_enu;
	  for(int ienu = 0; ienu < enu_spline->GetNbins(); ienu++){ // loop over true nu energy
#if USE_SPLINE_FD == USE_TSpline3_FD
		std::vector<TSpline3*> tmp_enu;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
		std::vector<TSpline3_red*> tmp_enu;
#endif
		std::vector<bool> tmp_flat_var1;
		std::vector<double> tmp_w_erec;
		for(int ierec = 0; ierec < var1_spline->GetNbins(); ierec++){ // loop over 1st variable
#if USE_SPLINE_FD == USE_TSpline3_FD
		  TSpline3 *tmp_erec=NULL;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
		  TSpline3_red *tmp_erec=NULL;
#endif
		  tmp_enu.push_back(tmp_erec);
		  tmp_flat_var1.push_back(false);//assume everything isn't flat intially
		  tmp_w_erec.push_back(1.0); // All weights can just be Set to 1
		} // end ierec loop
		tmp_mbin.push_back(tmp_enu);
		tmp_flat_enu.push_back(tmp_flat_var1);
		tmp_w_enu.push_back(tmp_w_erec);
	  } // end ienu loop               
	  tmp_tmp_imode.push_back(tmp_mbin);
	  tmp_flat_mode.push_back(tmp_flat_enu);
	  tmp_w_mode.push_back(tmp_w_enu);
	}// end of mode loop
	dev_1D_vec.push_back(tmp_tmp_imode);
	flat_vec.push_back(tmp_flat_mode);
	dev_1D_w.push_back(tmp_w_mode);
  }//end of syst loop

  for(int isyst=0 ; isyst < numSplineParams ; isyst++){
    syst* temp = new syst(SplineFileParsNames[isyst], &(dev_1D_vec.at(isyst)));
    systs.push_back(temp);
  }
 
  // Dummy spline: flat     
  TGraph *dummy_gr = new TGraph();
  dummy_gr->SetPoint(0,-99999999999,1);
  dummy_gr->SetPoint(1,0,1);
  dummy_gr->SetPoint(2,99999999999,1);

  /////////////////
  // Now load the splines from the spline file
  ////////////////

  TIter next(splinefile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TSpline3")) continue;

	char* splinename=(char*)key->GetName();
	//std::cout << "Spline is " << splinename << std::endl;
	char* syst;
	char* mode;
	int etruebin;
	int erecbin;

	char* tok= strtok(splinename,"_");//dev
	tok = strtok (NULL, "_");//syst
	syst = tok;

	int systnum=-1;
	for(unsigned isyst=0; isyst<systs.size(); isyst++){  // loop over systematics
	  if(strcmp(syst,systs.at(isyst)->name.c_str())==0){
		systnum=isyst;
		break;
	  }
	}

	//If the syst doesn't match any of the spline names then skip it
	//e.g. LowQ2suppression splines that we don't use in MaCh3
	if(systnum==-1){
	  continue;
	}	

	int modenum=-1;
	mode = strtok (NULL, "_");//mode
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
	  if(strcmp(mode,(UniqueModeFarSplineNames[imode]).c_str())==0){
		modenum=imode;
		break;
	  }
	}

	if(modenum==-1){
	  std::cout << "COULDN'T MATCH " << syst << std::endl;
	  std::cout << "No matching mode found for this spline... this shouldn't happen " << std::endl;
      throw;
	}

	tok = strtok (NULL, "_");//sp
	etruebin = atoi(strtok (NULL, "_"));//x
	erecbin = atoi(strtok (NULL, "_"));//y

	TSpline3 *h = (TSpline3*)key->ReadObj();
#if USE_SPLINE_FD == USE_TSpline3_FD
	TSpline3 *spl=(TSpline3*)h->Clone();
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	TSpline3_red *spl = new TSpline3_red(h);
#endif
	delete h;

	//loop over all the spline knots and check their value
	//if the value is 1 then Set the flat bool to false
	int n_knots = spl->GetNp();
	bool flat = true;
	for(int knot_i = 0 ; knot_i < n_knots ; knot_i++){
	  double x =-999;
	  double y = -999;
	  spl->GetKnot(knot_i, x, y);
	  if(x == -999 || y == -999){
		std::cerr << "Something has gone wrong... knot position is at -999" << std::endl;
		throw;
	  }
	  double eval = spl->Eval(x);
	  if(eval < 0.99999 || eval > 1.00001){flat = false; break;}
	}

	//If the spline is flat Set it to NULL and update the flat vector
	if(flat){
	  systs.at(systnum)->spline->at(modenum).at(etruebin).at(erecbin)=NULL; 
	  flat_vec.at(systnum).at(modenum).at(etruebin).at(erecbin) = flat;
	}
	else{
	  systs.at(systnum)->spline->at(modenum).at(etruebin).at(erecbin)=spl;
	}		
  }

  for(unsigned isyst=0; isyst<systs.size(); isyst++){  // loop over systematics
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
	  for(int ienu = 0; ienu < enu_spline->GetNbins(); ienu++){ // loop over true nu energy
		for(int ierec = 0; ierec < var1_spline->GetNbins(); ierec++){ // loop over 1st variable
		  bool flat_spline = flat_vec.at(isyst).at(imode).at(ienu).at(ierec); // check is the spline is flat
		  // if the spline is not flat and the spline is NULL then we have a problem!
		  if((systs.at(isyst)->spline->at(imode).at(ienu).at(ierec)==NULL) && (!flat_spline)){
			char sname[50];
			sprintf(sname,"dev_%s_%s_sp_%d_%d",systs.at(isyst)->name.c_str(),UniqueModeFarSplineNames[imode].c_str(),ienu,ierec);
			//Special case for params which apply to all modes i.e. I've Set mode = 12 in xsec cov 
			std::vector<int> modes = SplineModeVecs[isyst]; 
			for(unsigned spline_mode_i = 0 ; spline_mode_i < modes.size() ; spline_mode_i++){
			  if(modes[spline_mode_i] == imode){
				std::cerr << "[ERROR:] splinesDUNE::SetupSplines() - cannot FIND Erec SPLINE " << sname << std::endl;
				std::cerr << "[ERROR:] check that the spline name given in the xsec covariance matches that in the spline file" << std::endl;
			  }
			}//End of mode that splines apply to
		  }//End of if
		} // end ierec loop
	  } // end ienu loop               
	}//end of imode loop
  }//end of syst loop

  //splinefile->Close();                                       

  //We're now done with these structs so lets delete them
  for(unsigned int syst_i = 0 ; syst_i < systs.size() ; syst_i++){
    delete systs[syst_i];
  }

  return;
}


void splinesDUNE::SetupSplines(int opt_binning) // 2d version
{  
  BinningOpt = opt_binning;
  SetSplineBinning(BinningOpt); 

  int Nbins_1st_var=__BAD_SPLINE__, Nbins_2nd_var=__BAD_SPLINE__;

  Nbins_1st_var = var1_spline->GetNbins();
  Nbins_2nd_var = var2_spline->GetNbins();

  if (Nbins_1st_var == __BAD_SPLINE__ || Nbins_2nd_var == __BAD_SPLINE__) {
	std::cout << "Error: Nbins_1st_var or Nbins_2nd_var = __BAD_SPLINE__" << std::endl;
	std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
	exit(-1);
  }

  //DB Now use detid to determine number of spline systematics, names and corresponding modes
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);

  // ETA - need to think about how to do this configurably. If we store all the dev_blah_sp in one vector then can loop through giving the address of each?
  // Also need to pass in the name of the spline in the splinefile from the xsec xml to the covarianceXsec class
  std::vector<syst2D*> systs;

#if USE_SPLINE_FD == USE_TSpline3_red_FD
  std::cout << "###########################" << std::endl;
  std::cout << "USING TSPLINE3 RED !!!!!" << std::endl;
  std::cout << "###########################" << std::endl;
#endif

  std::vector<std::vector<std::vector<std::vector<std::vector<bool> > > > > flat_vec;
  for(int isyst=0; isyst< numSplineParams ; isyst++){  // loop over systematics
#if USE_SPLINE_FD == USE_TSpline3_FD
	std::vector<std::vector<std::vector<std::vector<TSpline3*> > > > tmp_tmp_imode;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	std::vector<std::vector<std::vector<std::vector<TSpline3_red*> > > > tmp_tmp_imode;
#endif
	std::vector<std::vector<std::vector<std::vector<bool> > > > tmp_flat_mode;
	//ETA adding this to store the weights for each spline eval
	std::vector<std::vector<std::vector<std::vector<double> > > > tmp_w_mode;
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
#if USE_SPLINE_FD == USE_TSpline3_FD
	  std::vector< std::vector<std::vector<TSpline3*> > > tmp_mbin;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	  std::vector< std::vector<std::vector<TSpline3_red*> > > tmp_mbin;
#endif
	  std::vector<std::vector<std::vector<bool> > > tmp_flat_enu;
	  std::vector<std::vector<std::vector<double> > > tmp_w_enu;
	  for(int ienu = 0; ienu < enu_spline->GetNbins(); ienu++){ // loop over true nu energy
#if USE_SPLINE_FD == USE_TSpline3_FD
		std::vector<std::vector<TSpline3*> > tmp_enu;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
		std::vector<std::vector<TSpline3_red*> > tmp_enu;
#endif
		std::vector<std::vector<bool> > tmp_flat_var1;
		std::vector<std::vector<double> > tmp_w_var1;
		for(int i1 = 0; i1 < Nbins_1st_var; i1++){ // loop over 1st variable
#if USE_SPLINE_FD == USE_TSpline3_FD
		  std::vector<TSpline3*> tmp_var1;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
		  std::vector<TSpline3_red*> tmp_var1;
#endif
		  std::vector<bool> tmp_flat_var2;
		  std::vector<double> tmp_w_var2;
		  for (int i2 = 0; i2 < Nbins_2nd_var; i2++){ // loop over 2nd variable
#if USE_SPLINE_FD == USE_TSpline3_FD
			TSpline3 *tmp_var2=NULL;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
			TSpline3_red *tmp_var2=NULL;
#endif
			tmp_var1.push_back(tmp_var2);
			tmp_flat_var2.push_back(false);//assume everything is flat initally
			tmp_w_var2.push_back(1.0); // All weights can just be Set to 1
		  } // end i2 loop
		  tmp_enu.push_back(tmp_var1);
		  tmp_flat_var1.push_back(tmp_flat_var2);
		  tmp_w_var1.push_back(tmp_w_var2);
		} // end i1 loop
		tmp_mbin.push_back(tmp_enu);
		tmp_flat_enu.push_back(tmp_flat_var1);	
		tmp_w_enu.push_back(tmp_w_var1);
	  } // end ienu loop               
	  tmp_tmp_imode.push_back(tmp_mbin);
	  tmp_flat_mode.push_back(tmp_flat_enu);
	  tmp_w_mode.push_back(tmp_w_enu);
        }// end of mode loop
	flat_vec.push_back(tmp_flat_mode);
	dev_2D_vec.push_back(tmp_tmp_imode);
	dev_2D_w.push_back(tmp_w_mode);
      }// end of syst loop

  for(int isyst=0 ; isyst < numSplineParams ; isyst++){
     	syst2D* temp = new syst2D(SplineFileParsNames[isyst], &(dev_2D_vec.at(isyst)));
       	systs.push_back(temp);
       }

  // Dummy splines: flat               
  TGraph *dummy_gr = new TGraph();
  dummy_gr->SetPoint(0,-99999999999,1);
  dummy_gr->SetPoint(1,0,1);
  dummy_gr->SetPoint(2,99999999999,1);

  std::cout<<"dummy gr made"<<std::endl;
  /////////////////
  // Now load the splines from the spline file
  ////////////////
    
  //get all spline objects from file
  TIter next(splinefile->GetListOfKeys());
  std::cout<<"GetListOfKeys"<<std::endl;
  TKey *key;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TSpline3")) continue;

    //std::cout<< "Spline is " << key->GetName()<<std::endl; 

    char* splinename=(char*)key->GetName();
    char* syst;
    char* tok= strtok(splinename,"_");//dev
    tok = strtok (NULL, "_");//syst
	syst=tok;
    

    int systnum=-1;
    for(unsigned isyst=0; isyst<systs.size(); isyst++){  // loop over systematics
      //std::cout<<"looping over systematics, isyst: "<<isyst<<std::endl;
      if(strcmp(syst,systs.at(isyst)->name.c_str())==0){
        systnum=isyst;
	break;
      }
    }
    
    //If the syst doesn't match any of the spline names then skip it
    //e.g. LowQ2suppression splines that we don't use in MaCh3
    if(systnum==-1){
      continue;
    }
    
    int modenum=-1;
    char* mode = strtok (NULL, "_");//mode
    for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
      if(strcmp(mode,(UniqueModeFarSplineNames[imode]).c_str())==0){
	modenum=imode;
	break;
      }
    }
    
    if(modenum==-1){
      std::cout << "COULDN'T MATCH " << syst << std::endl;
      std::cout << "No matching mode found for this spline... this shouldn't happen " << std::endl;
      throw;
    }
    
    tok = strtok (NULL, "_");//sp
    int etruebin = atoi(strtok (NULL, "_"));//x
    int var1bin = atoi(strtok (NULL, "_"));//y
    int var2bin = atoi(strtok (NULL, "_"));//z
 
    TSpline3 *h = (TSpline3*)key->ReadObj();

#if USE_SPLINE_FD == USE_TSpline3_FD
    TSpline3 *spl=(TSpline3*)h->Clone();
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	TSpline3_red *spl = new TSpline3_red(h);
#endif

	delete h;
	//std::cout << "address is " << spl << std::endl;
	
	//loop over all the spline knots and check their value
	//if the value is 1 then Set the flat bool to false
	int n_knots = spl->GetNp();
	bool flat = true;
	for(int knot_i = 0 ; knot_i < n_knots ; knot_i++){
          //std::cout<<"looping over knots"<<std::endl;

	  double x =-999;
	  double y = -999;
	  spl->GetKnot(knot_i, x, y);
	  if(x == -999 || y == -999){
		std::cerr << "Something has gone wrong... knot position is at -999" << std::endl;
		throw;
	  }
	  double eval = spl->Eval(x);
	  if(eval < 0.99999 || eval > 1.00001){flat = false; break;}
	}

	//If the spline is flat Set it to NULL and update the flat vector
	if(flat){
	  systs.at(systnum)->spline->at(modenum).at(etruebin).at(var1bin).at(var2bin)=NULL;
	  flat_vec.at(systnum).at(modenum).at(etruebin).at(var1bin).at(var2bin) = true;
	}
	else{
	  systs.at(systnum)->spline->at(modenum).at(etruebin).at(var1bin).at(var2bin)=spl;
	} 

  } 

  for(unsigned isyst=0; isyst<systs.size(); isyst++){  // loop over systematics
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
	  for(int ienu = 0; ienu < enu_spline->GetNbins(); ienu++){ // loop over true nu energy
		for(int i1 = 0; i1 < Nbins_1st_var; i1++){ // loop over 1st variable
		  for (int i2 = 0; i2 < Nbins_2nd_var; i2++){ // loop over 2nd variable
			bool flat_spline = flat_vec.at(isyst).at(imode).at(ienu).at(i1).at(i2);//check to see if the spline is flat
		  // if the spline is not flat and the spline is NULL then we have a problem!
			if((systs.at(isyst)->spline->at(imode).at(ienu).at(i1).at(i2)==NULL) && (!flat_spline)){
			  char sname[50];
			  sprintf(sname,"dev_%s_%s_sp_%d_%d_%d",systs.at(isyst)->name.c_str(),UniqueModeFarSplineNames[imode].c_str(),ienu,i1,i2);
			  //Special case for params which apply to all modes i.e. I've Set mode = 12 in xsec cov 
			  std::vector<int> modes = SplineModeVecs[isyst];
			  for(unsigned spline_mode_i = 0 ; spline_mode_i < modes.size() ; spline_mode_i++){
				if(modes[spline_mode_i] == imode){
				  std::cerr << "[ERROR:] splinesDUNE::SetupSplines() - cannot FIND Erec SPLINE " << sname << std::endl;
				  std::cerr << "[ERROR:] check that the spline name given in the xsec covariance matches that in the spline file" << std::endl;
				}
			  }//End of mode that splines apply to
			}
		  } // end i2 loop
		} // end i1 loop
	  } // end ienu loop               
	}
  }

//  splinefile->Close();      


  //We're now done with these structs so lets delete them
  for(unsigned int syst_i = 0 ; syst_i < systs.size() ; syst_i++){
    //std::cout<<"deleting syst_i: "<<syst_i<<std::endl;

    delete systs[syst_i];
  }
  std::cout<<"deleted all structs in systs vector "<<std::endl;

  return;
}

//TODO (ETA) - need to pass number of interaction modes and unique spline modes to spline object
//the mode inofrmation etc. will be defined for each experiment
void splinesDUNE::FindUniqueModes() {
  
  nUniqueModes = 0;

  std::cout << "in FindUniqueModes " << std::endl;

  for (int iMode=0;iMode<kMaCh3_nModes;iMode++) {
    if (MaCh3Mode_to_SplineMode(iMode)==iMode) {
      nUniqueModes += 1;
    }
  }

  int Counter = 0;
  UniqueModeFarSplineNames.resize(nUniqueModes);

  for (int iMode=0;iMode<kMaCh3_nModes;iMode++) {
    if (MaCh3Mode_to_SplineMode(iMode)==iMode) {
      UniqueModeFarSplineNames[Counter] = MaCh3mode_ToDUNEString((MaCh3_Mode)iMode);
    } else {
      DuplicatedFDModes.push_back(iMode);
    }
    Counter += 1;
  }

  MaCh3Mode_SplineMode_Map.resize(kMaCh3_nModes);
  for (int iMode=0;iMode<kMaCh3_nModes;iMode++) {
    MaCh3Mode_SplineMode_Map[iMode] = MaCh3Mode_to_SplineMode(iMode);
  }
   

/* for (int i = 0; i < 1;i++) {
   MaCh3Mode_SplineMode_Map.push_back(i); }
   nUniqueModes = 1;
   UniqueModeFarSplineNames.push_back("ccqe"); 
   std::cout << "Temp implementation" << std::endl;
   */

  return;
}
