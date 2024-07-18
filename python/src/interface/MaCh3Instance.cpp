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
    std::vector<double> oscpars = fit_manager->raw()["General"]["OscParams"].as<std::vector<double>>();

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
    std::vector<std::string> fd_sample_configs = fit_manager->raw()["General"]["FDSamples"].as<std::vector<std::string>>();
    std::vector<std::string> nd_sample_configs = fit_manager->raw()["General"]["NDSamples"].as<std::vector<std::string>>();

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


void MaCh3Instance::set_parameter_values(std::vector<double> new_pars){
    // This is HACKED together
    if((int)new_pars.size()<osc->getNpars()+xsec->getNpars()){
        std::cerr<<"ERROR::Trying to set wrong number of values"<<std::endl;
        std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
        throw;
    }

    for(int i=0; i<(int)new_pars.size(); i++){
        if(i<xsec->getNpars()){
            xsec->setParProp(i, new_pars[i]);
        }
        else{
            osc->setParProp(i-xsec->getNpars(), new_pars[i]);
        }

    }
}

double MaCh3Instance::get_likelihood(){
    double llh = 0;

    llh += osc->GetLikelihood();
    llh += xsec->GetLikelihood();

    // Reject based on boundary conditions
    if(llh>=__LARGE_LOGL__){
        return llh * sample_vector.size(); // scale by NSamples
    }

    // reweight samples
    for(const auto sample : sample_vector){
        sample->reweight(osc->getPropPars());
        llh+=sample->GetLikelihood();
    }

    return llh;

}

std::vector<double> MaCh3Instance::get_parameter_values(){
    //hacked in method for getting parameters
    std::vector<double> xsec_pars = xsec->getParameters();
    std::vector<double> osc_pars = osc->getParameters();

    xsec_pars.insert(xsec_pars.end(), osc_pars.begin(), osc_pars.end());
    return xsec_pars;
}

double MaCh3Instance::propose_step(std::vector<double> new_step){
    set_parameter_values(new_step);
    return get_likelihood();
}


// ################################################################################################
// PyBind wrapping
PYBIND11_MODULE(_pyMaCh3, m){
    py::class_<MaCh3Instance>(m, "MaCh3Instance")
        .def(py::init<std::string>())
        .def("get_parameter_values", &MaCh3Instance::get_parameter_values)
        .def("propose_step", &MaCh3Instance::propose_step);
}
