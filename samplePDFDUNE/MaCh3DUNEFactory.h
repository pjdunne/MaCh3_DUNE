#pragma once

// Include the input manager
#include "manager/manager.h"

// Include the samplePDFs
#include "samplePDFDUNE/samplePDFDUNEBeamFDBase.h"

void MakeMaCh3DuneBeamInstance(manager *fitMan, std::vector<samplePDFDUNEBeamFDBase*> &sample_vec,  covarianceXsec *&xsec, covarianceOsc *&osc);
