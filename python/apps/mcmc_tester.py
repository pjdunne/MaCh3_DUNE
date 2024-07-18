'''
MCMC interface test with MaCh3!
'''
from _pyMaCh3 import MaCh3Instance
import argparse
import numpy as np

import tensorflow as tf
import tensorflow_probability as tfp



if __name__=="__main__":
    parser = argparse.ArgumentParser(usage="python demo_app -c <config_name>.yaml")
    parser.add_argument("-c", "--config", help="YAML config file")
    
    args = parser.parse_args()

    mach3 = MaCh3Instance(args.config)
    initial_values = np.transpose(mach3.get_parameter_values())
    
    
    def custom_get_likelihood(x):
        x_list = x.numpy().tolist()
        return mach3.propose_step(x_list)

    
    num_results = int(10e3)
    num_burnin_steps = int(1e3)
    adaptive_hmc = tfp.mcmc.SimpleStepSizeAdaptation(
    tfp.mcmc.HamiltonianMonteCarlo(
        target_log_prob_fn=custom_get_likelihood,
        num_leapfrog_steps=3,
        step_size=1.),
        num_adaptation_steps=int(num_burnin_steps * 0.8)
    )

    @tf.function
    def run_chain():
        # Run the chain (with burn-in).
        samples, is_accepted = tfp.mcmc.sample_chain(
            num_results=num_results,
            num_burnin_steps=num_burnin_steps,
            current_state=initial_values,
            kernel=adaptive_hmc,
            trace_fn=lambda _, pkr: pkr.inner_results.is_accepted)

        sample_mean = tf.reduce_mean(samples)
        sample_stddev = tf.math.reduce_std(samples)
        is_accepted = tf.reduce_mean(tf.cast(is_accepted, dtype=tf.float32))
        return sample_mean, sample_stddev, is_accepted


    sample_mean, sample_stddev, is_accepted = run_chain()

    print('mean:{:.4f}  stddev:{:.4f}  acceptance:{:.4f}'.format(
        sample_mean.numpy(), sample_stddev.numpy(), is_accepted.numpy()))
