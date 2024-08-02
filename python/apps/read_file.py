'''
HW : Simple script for reading and plotting the diagnostics of plots created by emcee
'''

import emcee
import matplotlib.pyplot as plt
import argparse
import numpy as np
import corner

if __name__=="__main__":
    parser = argparse.ArgumentParser(usage="python read_file -i <input_file>.h5 Optional -o <output_file>.pdf")
    parser.add_argument("-i", "--input_file", help="H5 input file")
    parser.add_argument("-o", "--output_file", help="PDF output file")

    args = parser.parse_args()

    reader = emcee.backends.HDFBackend(args.input_file)
    
    tau = reader.get_autocorr_time()
    burnin = int(2 * np.max(tau))
    thin = int(0.5 * np.min(tau))
    samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
    log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
    log_prior_samples = reader.get_blobs(discard=burnin, flat=True, thin=thin)

    print(f"burn-in: {burnin}")
    print(f"thin: {thin}")
    print(f"flat chain shape: {samples.shape}")
    print(f"flat log prob shape: {log_prob_samples.shape}")
    print(f"flat log prior shape: {log_prob_samples.shape}")
    
    