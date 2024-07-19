#pragma once

#include "pyMaCh3.h"
// C++ includes
#include <string>
#include <vector>
#include <type_traits>

// MaCh3 Includes
#include "manager.h"
#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include "samplePDFDUNE/samplePDFDUNEBaseND.h"

// pybind
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

// pyMaCh3
#include "pyMaCh3.h"


namespace py = pybind11;

class MaCh3SampleHandler{
 public:
    /// @brief Sample handler constructor
    MaCh3SampleHandler();

    /// @brief gets number of events in each data bin
    std::vector<double> get_data_bin_values();

    /// @brief gets number of events in each MC bin [uses same binning loop as data]
    std::vector<double> get_mc_bin_values();


 protected:
    
    /// @brief sets up samples from vector of names taken from YAML config
    template<class T>
    void setup_samples(std::vector<std::string> sample_names, double pot){
        
        for(std::string sample : sample_names){
            static_assert(std::is_base_of<samplePDFFDBase, T>::value, "If Error: samplePDFFDBase is not base of T");
            
            T* sample_pdf = new T(pot, sample.c_str(), xsec);

            sample_pdf->useNonDoubledAngles(true);
    
            sample_pdf->SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
            
            // Add Data
            sample_pdf->reweight(osc->getPropPars());

            // Get name
            TString samp_name(sample_pdf->GetSampleName().c_str());

            switch(sample_pdf->GetBinningOpt()){
                case 1:
                    sample_pdf->addData((TH1D*)sample_pdf->get1DHist()->Clone(samp_name+"_asimov"));
                    break;
                case 2:
                    sample_pdf->addData((TH2D*)sample_pdf->get2DHist()->Clone(samp_name+"_asimov"));
                    break;
                default:
                    std::cerr<<"ERROR::Unknown sample binning opt, please pick 1 (1D) or 2 (2D)"<<std::endl;
                    std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
                    throw;
            }

            sample_vector.push_back(sample_pdf);
        }
    }

    std::vector<samplePDFFDBase*> sample_vector;

};