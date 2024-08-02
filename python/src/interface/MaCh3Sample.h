#pragma once

// MaCh3 includes
#include "samplePDF/samplePDFBase.h"
#include "MaCh3Object.h"

// C++ Includes
#include <iostream>
#include <vector>
#include <string>

class MaCh3Systematic : public MaCh3Object<samplePDFBase>{
   public:
    /// @brief MaCh3Systematic constructor
    MaCh3Systematic(std::string config_opt);

    /// @brief destructor
    virtual ~MaCh3Systematic(){}

    /// @brief  Likelihood getter
    double get_likelihood(){
        return _stored_object->GetLikelihood();
    }

    /// @brief wrapper around the reweight method
    void reweight(std::vector<double> osc_parameter_values){
        // Reweight takes double*
        _stored_object->reweight(osc_parameter_values.data());
    }

   protected:

    /// @brief Constructs dummy samplePDFBase, since this is never actually likely to be used I'm giving it a dummy arg!
    std::unique_ptr<samplePDFBase> construct_default() override{
        return std::unique_ptr<samplePDFBase>(samplePDFBase(0));
    }


};