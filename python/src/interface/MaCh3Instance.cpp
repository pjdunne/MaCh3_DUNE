#include "MaCh3Instance.h"

// Public methods
MaCh3Instance::MaCh3Instance(std::string yaml_config){
    // Open YAML config with manager wrapper
    manager *fit_manager = new manager(yaml_config);

    // Lets grab basic information from config
    double fd_pot = fit_manager->raw()["General"]["FDPOT"].as<double>(); 
    double nd_pot = fit_manager->raw()["General"]["NDPOT"].as<double>(); 

    // Load Xsec Matrix
    std::string  xsec_matrix_file = fit_manager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>(); 
    std::string  xsec_matrix_name = fit_manager->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
    
    xsec = new covarianceXsec(xsec_matrix_name.c_str(), xsec_matrix_file.c_str());
    xsec->setParameters();


    std::string  osc_matrix_file = fit_manager->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
    std::string  osc_matrix_name = fit_manager->raw()["General"]["Systematics"]["OscCovName"].as<std::string>();
  
    osc = new covarianceOsc(osc_matrix_name.c_str(), osc_matrix_file.c_str());

    if (!(oscpars.size()==6 || oscpars.size()==7))
    {
        std::cout<<"Input osc pars not of right size, there should be six entries (or seven if setting beta)"<<std::endl;
        std::cout<<"oscpars.size() = " << oscpars.size() << std::endl;
        exit(1);
    }

    std::cout<<"Using these oscillation parameters: ";
    for(unsigned ipar=0;ipar<oscpars.size();ipar++){
        std::cout<<" "<<oscpars.at(ipar);
    }
    std::cout << std::endl;


    osc->setFlipDeltaM23(true);

    // Use prior for 12 parameters only
    //osc->setEvalLikelihood(0,false);
    osc->setEvalLikelihood(1,false);
    osc->setEvalLikelihood(2,false);
    //osc->setEvalLikelihood(3,false);
    osc->setEvalLikelihood(4,false);
    osc->setEvalLikelihood(5,false);

    osc->setParameters(oscpars);
    osc->acceptStep();


    // Load in samples from config
    std::vector<std::string> fd_sample_configs = FitManager->raw()["General"]["FDSamples"].as<std::vector<std::string>>();
    std::vector<std::string> nd_sample_configs = FitManager->raw()["General"]["NDSamples"].as<std::vector<std::string>>();

    // This will initialise all samples
    setup_samples<samplePDFDUNEBase>(fd_sample_configs, fd_pot);
    setup_samples<samplePDFDUNEBaseND>(fd_sample_configs, nd_pot);

    // Check this has done something!
    if(sample_vector.size()==0){
        std::cerr<<"ERROR::No samples have been initialised!"<<std::endl;
        std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
        throw;
    }



    xsec->proposeStep();
    osc->proposeStep();
    
    
    // set up parameter indices :: NOTE, this is hacked in and fixed in later versions of core

    std::cout<<"All samples+systematics intialised, ready to roll!"<<std::endl;
}

void MaCh3Interface::set_parameter_values(std::vector<double> new_pars){
    // This is HACKED together
    if((int)new_pars.size()<osc->getNpars()+xsec->getNpars()){
        std::cerr<<"ERROR::Trying to set wrong number of values"<<std::endl;
        std::cerr<<__FILE__<<":"<<__LINE<<std::endl;
        throw;
    }

    for(unsigned i=0; i<new_pars.size(); i++){
        if(i<xsec.size()){
            xsec->setParProp(i, new_pars[i]);
        }
        else{
            osc->setParProp(i-xsec.size(), new_pars[i]);
        }

    }
}

double MaCh3Interface::get_likelihood(){
    double llh = 0;
    reject = false;

    llh += osc->GetLikelihood();
    llh += xsec->GetLikelihood();

    // Reject based on boundary conditions
    if(llh>=__LARGE_LOGL__){
        return llh * sample_vector.size(); // scale by NSamples
    }

    // reweight samples
    for(const auto sample : sample_vector){
        sample->reweight(ocs->getPropPars());
        llh+=sample->GetLikelihood();
    }

    return llh;

}

double propose_step(std::vector<double> new_step){
    set_parameter_values(new_step);
    return get_likelihood();
}


// ################################################################################################
// Protected
template<typename T> 
void MaCh3Interface::setup_samples(std::vector<std::string> sample_names, double pot){
    for(std::string sample : sample_names){
        static_assert(is_base<samplePDFFDBase,T>::value, "If Error: samplePDFFDBase is not base of T");
        
        T* sample_pdf = new T(pot, sample.c_str(), xsec);

        sample_pdf->useNonDoubleAngles(true);
        sample_pdf->SetupOscCalc(osc->GetPathLength(), osc->GetDensity());

        // Add Data
        sample_pdf->reweight(osc->getPropPars());

        // Get name
        TString samp_name(sample_pdf->GetSampleName().c_str());

        switch(sample_pdf->GetBinningOpt()){
            case 1:
                sample_pdf->addData((TH1D*)SamplePDFs[sample_i]->get1DHist()->Clone(samp_name+"_asimov"));
                break;
            case 2:
                sample_pdf->addData((TH3D*)SamplePDFs[sample_i]->get2DHist()->Clone(samp_name+"_asimov"));
                break;
            default:
                std::cerr<<"ERROR::Unknown sample binning opt, please pick 1 (1D) or 2 (2D)"<<std::endl;
                std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
                throw;
        }

        sample_vector.push_back(sample_pdf);
    }
}


