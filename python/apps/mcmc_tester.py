'''
MCMC interface test with MaCh3!
'''
from pyMaCh3 import MaCh3Instance
import argparse
import numpy as np
from typing import Iterable
import emcee
from functools import partial
import matplotlib.pyplot as plt
from jax import value_and_grad


#HACK: Wrapper around propose step, should really be done at the pybind level!
def propose_step(mach3_instance: MaCh3Instance, input_iterable: Iterable):
    return -1*mach3_instance.propose_step(list(input_iterable))

if __name__=="__main__":

    parser = argparse.ArgumentParser(usage="python demo_app -c <config_name>.yaml Optional :[ -w <number of walkers> -n <number of steps]")
    parser.add_argument("-c", "--config", help="YAML config file")
    parser.add_argument("-w", "--n_walkers", default=20, type=int, help="Number of walkers")
    parser.add_argument("-n", "--n_steps", default=100, type=int, help="Number of Steps")
    
    args = parser.parse_args()

    mach3 = MaCh3Instance(args.config)

    initial_values = mach3.get_parameter_values()

    likelihood_func = partial(propose_step, mach3)

    ndim = len(initial_values)
    
    val, gradient = value_and_grad(likelihood_func(initial_values+0.0001*np.random.rand(ndim)))
    
    print(f"Value: {val}   | Grad: {gradient}")
    
    # print(f"Dimension is {ndim}, using {args.n_walkers} walkers")
    # p0 = [initial_values+0.0001*np.random.rand(ndim) for _ in range(args.n_walkers)]
    
    # filename = "tutorial.h5"
    # backend = emcee.backends.HDFBackend(filename)
    # backend.reset(args.n_walkers, ndim)

    
    
    # my_sampler = emcee.EnsembleSampler(args.n_walkers, ndim, likelihood_func, live_dangerously=True, backend=backend, threads=8)
    # my_sampler.run_mcmc(p0, args.n_steps, skip_initial_state_check=True, progress=True)
    # print(f"Accepted {my_sampler.acceptance_fraction[-1]}% of steps")
    
    
    # samples = my_sampler.get_chain(flat=True, discard = int(0.01*args.n_steps))
    
    
    
    # plt.hist(samples[:, 0], 100, color="k", histtype="step")
    # plt.xlabel(r"$\theta_1$")
    # plt.ylabel(r"$p(\theta_1)$")
    # plt.gca().set_yticks([]);

    # plt.savefig("test_fig_2.png")
    # plt.close()