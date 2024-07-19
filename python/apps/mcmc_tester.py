'''
MCMC interface test with MaCh3!
'''
from pyMaCh3 import MaCh3Instance
import argparse
import numpy as np
from typing import Iterable
import emcee
from functools import partial


#HACK: Wrapper around propose step, should really be done at the pybind level!
def propose_step(mach3_instance: MaCh3Instance, input_iterable: Iterable):
    return mach3_instance.propose_step(list(input_iterable))

if __name__=="__main__":

    parser = argparse.ArgumentParser(usage="python demo_app -c <config_name>.yaml Optional :[ -w <number of walkers>]")
    parser.add_argument("-c", "--config", help="YAML config file")
    parser.add_argument("-w", "--n_walkers", default=20, type=int, help="Number of walkers")
    
    args = parser.parse_args()

    mach3 = MaCh3Instance(args.config)

    initial_values = mach3.get_parameter_values()

    likelihood_func = partial(propose_step, mach3)
    
    ndim = len(initial_values)
    print(f"Dimension is {ndim}, using {args.n_walkers}")
    p0 = [initial_values+0.01*np.random.rand(ndim) for _ in range(args.n_walkers)]
    
    
    my_sampler = emcee.EnsembleSampler(args.n_walkers, ndim, likelihood_func, live_dangerously=True)
    my_sampler.run_mcmc(p0, 10000, skip_initial_state_check=True, progress=True)