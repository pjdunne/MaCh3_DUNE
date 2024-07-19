'''
MCMC interface test with MaCh3!
'''
from pyMaCh3 import MaCh3Instance
import argparse
import numpy as np

import emcee

if __name__=="__main__":

    parser = argparse.ArgumentParser(usage="python demo_app -c <config_name>.yaml Optional :[ -w <number of walkers>]")
    parser.add_argument("-c", "--config", help="YAML config file")
    parser.add_argument("-w", "--n_walkers", default=20, help="Number of walkers")
    
    args = parser.parse_args()

    mach3 = MaCh3Instance(args.config)

    initial_values = np.array(mach3.get_parameter_values())
    
    
    args.n_walkers = 3
    ndim = len(initial_values)
    p0 = [initial_values+0.01*np.random.rand(ndim) for _ in range(args.n_walkers)]
    
    print(p0)
    sampler = emcee.EnsembleSampler(args.n_walkers, ndim, mach3.propose_step, skip_initial_state_check=True)
    
    sampler.run_mcmc(p0, 100)