#pragma once

/*
Oscillation Systematic plugin
*/

// MaCh3 includes
#include "../src/interface/MaCh3Plugin.h"
#include "../src/interface/MaCh3Systematic.h"
#include "manager/manager.h"

#include "covariance/covarianceOsc.h"


#include <memory>
#include <vector>#

class OscSystPlugin : public MaCh3Plugin<covarianceOsc>{
    OscSystPlugin() : MaCh3Plugin() {
        override_default_constructor = true;
    }

    std::unique_ptr<covarianceOsc> construct_obj(std::string matrix_name, std::string file_name){
        auto Osc = std::make_unique<covarianceOsc>(matrix_name, file_name);
        return Osc;
    }


    void process(MaCh3Systematic * mach3_object){
        mach3_object->get_stored_object()->setFlipDeltaM23(true);
        mach3_object->get_stored_object()->setEvalLikelihood(1,false);
        mach3_object->get_stored_object()->setEvalLikelihood(2,false);
        mach3_object->get_stored_object()->setEvalLikelihood(4,false);
        mach3_object->get_stored_object()->setEvalLikelihood(5,false);

        std::vector<double> oscpars = mach3_object->get_manager()->raw()["General"]["OscParams"].as<std::vector<double>>();
        
        mach3_object->get_stored_object()->setParameters(oscpars);
        mach3_object->get_stored_object()->acceptStep();

    }
};

