# pyMaCh3
NOTE All work here is preliminary and may break!

## Introduction
The following code provides an API for interfacing between python scripts and MaCh3.

## Functionality
Currently the API can do 3 things 
- Create a MaCh3 instance (load inputs etc.) from a standard MaCh3 YAML config
- Get the current parameter values from MaCh3
- Get the likelihood for a given set of parameter values
 
## How Do I run this?
The python libraries are stored in `build/python/pyMaCh3` and can be accessed like any usual python library. Examples of this can be found in `python/apps/demo_app.py` and `python/apps/mcmc_tester.py`. Running is very simple

Source necessary scripts:
```bash
source /path/to/this_root.sh # Needed for correct pathing!
source build/bin/setup.MaCh3DUNE.sh
```
I'd recommend using a virtual enivironment
```bash
virtualenv --no-download pyMaCh3_ENV
source pyMaCh3_ENV/bin/activate
```
The python depencies for running the 2 scripts can then be installed with
```bash
pip install -r python/requirements.txt
```
Running the scripts is then simple. In your top MaCh3_DUNE directory run
```
python build/python/demo_app.py  -c /path/to/yaml 
```
This is a simple app that just varies the value of one parameter by a small amount and outputs the likelihood. 

To run a a Markov Chain you can try
```bash
python build/python/mcmc_tester.py -c /path/to/yaml -n <number of steps> -w <number of walkers>
```
This runs an Emcee [https://emcee.readthedocs.io/en/stable/] chain using the MaCh3 likelihood handler

## Guide for Devs
<tbf>