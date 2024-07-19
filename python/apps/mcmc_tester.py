'''
MCMC interface test with MaCh3!
'''
from _pyMaCh3 import MaCh3Instance
import argparse
import numpy as np

import emcee

if __name__=="__main__":

    parser = argparse.ArgumentParser(usage="python demo_app -c <config_name>.yaml")
    parser.add_argument("-c", "--config", help="YAML config file")
    
    args = parser.parse_args()

    mach3 = MaCh3Instance(args.config)

    initial_values = np.array(mach3.get_parameter_values())
    
    
    nwalkers = 3
    ndim = len(initial_values)
    p0 = [[initial_values+0.01*np.random.rand(ndim)] for _ in range(nwalkers)]
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, mach3.propose_step)
    
    sampler.run_mcmc(p0, 100)