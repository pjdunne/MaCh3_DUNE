from pyMaCh3 import pyMaCh3Instance
import argparse

# HW : Really really simple python app

if __name__=="__main__":
    parser = argparse.ArgumentParser(usage="python demo_app -c <config_name>.yaml")
    parser.add_argument("-c", "--config", help="YAML config file")
    
    args = parser.parse_args()
    
    mach3 = pyMaCh3Instance(args.config)
    
    print(f"Parameter values are : {mach3.get_parameter_values}")
    
    parameter_values = mach3.get_parameter_values()
    parameter_values[0] += 1
    print(f"Likelihood after shift : {mach3.propose_step(parameter_values)}")
    

    