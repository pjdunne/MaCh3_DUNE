#include "TString.h"
#include "OscClass/OscClass_CUDAProb3.h"

int main(int argc, char const *argv[])
{
    // TString earthFile("/home/pgranger/atmospherics/MaCh3_DUNE/build/_deps/cudaprob3-src/models/PREM_4layer.dat");
    std::string config("/home/pgranger/atmospherics/debug/MaCh3_DUNE/configs/OscillatorObj.yaml");
    Oscillator osc(config);
    double oscpars[6] = {0.307,0.528,0.0218,7.53e-5, 2.509e-3,-1.601};
    double prodH = 20;
    std::cout << "Going to fill oscillogram" << std::endl;
    osc.SetFillHistograms();
    // osc.SetFillPrimaryOscillogram();
    osc.FillOscillogram(oscpars, prodH);
    osc.SaveOscillogramsToFile("test.root");
    /* code */
    return 0;
}
