import subprocess
import os
from ruamel.yaml import YAML
import argparse

def modify_config(config_file_path, searchtype, searchstring, newstring, KinematicStr=True, SamplePDF=True):
    yaml=YAML()
    with open(config_file_path) as f:
        conf=yaml.load(f)
#        yaml.default_flow_style=False
#        yaml_lines = yaml.dump(conf)
#    print(conf)

#    print(len(yaml_lines))
#    print(yaml_lines[2])
# Find the index of the line after the search string
    found = False
    index_to_modify =0
    if SamplePDF:
        if KinematicStr:
            for obj in conf[searchtype]:
#           print(obj)
                if obj['KinematicStr'] == searchstring:
 #                  print(len(obj))
                    obj['Bounds'] = newstring
 #                  print(newstring)
                    break
        if args.withecalcontainment:
            conf['SampleBools']['ecal_containment'] = 'yes'
        elif not args.withecalcontainment:
            conf['SampleBools']['ecal_containment'] = 'no'
        if not args.varyfdv:
            conf['SampleCuts']['TPCFidRadius'] = args.tpcfidrad
            conf['SampleCuts']['TPCInstrumentedRadius'] = args.tpcinstrumentedrad
            conf['SampleCuts']['ECALInnerRadius'] = ecalinnerrad
            conf['SampleCuts']['ECALOuterRadius'] = ecalouterrad
        elif args.varyfdv:
            if searchstring == 'TPCFidRad':
                conf['SampleCuts']['TPCFidRadius'] = newstring 
        if not args.varyBfield:
            conf['SampleCuts']['B_field']=args.bfield
        elif args.varyBfield:
            if searchstring == 'B_field':
                conf['SampleCuts']['B_field']=newstring
    elif not SamplePDF:
        conf['General']['Output']['FileName'] = newstring
        print(conf)

    with open(config_file_path, 'w') as f:
        yaml.dump(conf, f)
#    for i, line in enumerate(yaml_lines):
#        if searchstring in line:
#            found = True
#            if nextline:
#                index_to_modify = i + 1
#                break
#            else:
#                index_to_modify = i
#                break
#    print(index_to_modify)
#    if found:
#        yaml_lines[index_to_modify] = newstring
#        print(yaml_lines[1])
#        # Write the modified lines back to the YAML file
#        with open(config_file_path, 'w') as f:
#            yaml.dump(yaml_lines)


#    else:
#        print(f"Search string '{searchstring}' not found in the YAML file.")

        
def main():
    if args.withecalcontainment:
        ecalcontainment = "withecalcontainment"
    else:
        ecalcontainment = "noecalcontainment"
    if args.varyenergy:
        enu = [0.00, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 4.00, 4.50, 5.00]
#        enu=[1.50, 1.75]
        accepted_arr=[0.0, 1.0, 2.0]
        print(len(accepted_arr))
        for i in range(0, len(enu)-1):
            config_samplepdf = 'configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml'
            kinematicstr = 'TrueNeutrinoEnergy'
#            newvalues = '[{boundlow:.2f}, {boundup:.2f}]'.format(boundlow=enu[i], boundup=enu[i+1])
            newvalues = [enu[i], enu[i+1]]
            modify_config(config_samplepdf, 'SelectionCuts', kinematicstr, newvalues, True)
            modify_config(config_samplepdf, 'SelectionCuts', kinematicstr, newvalues, True)
            for j in range(0, len(accepted_arr)-1):
                print(accepted_arr)
                print(j)
                modify_config(config_samplepdf, 'SelectionCuts', 'IsAccepted', [accepted_arr[j], 2.0], True) 
                config_selec = 'configs/Selections_NDGAr_acceptancecorrection_2dhist.yaml'
                nameout = 'Output: '
                if j == 1:
                    accepted = "Accepted"
                else:
                    accepted = "All"
                filename = 'outputs/01_12_2024/Selections_CC_AcceptanceCorrection_2dhist_{accepted}Events_2point5mmpointres_fidrad{fidrad}_Bfield{bfield}_Enubin{binnum}_{ecal_containment}'.format(accepted=accepted, fidrad=str(args.tpcfidrad).replace(".", "_"), bfield = str(args.bfield).replace(".", "_"), binnum=i, ecal_containment=ecalcontainment)
                newoutput = '{filename}.root'.format(filename=filename)
                modify_config(config_selec, 'Output', nameout, newoutput, False, False)
                command = 'Selections configs/Selections_NDGAr_acceptancecorrection_2dhist.yaml &> {filename}.txt'.format(filename=filename)
                print(command)
                p = subprocess.Popen(command, shell=True)
                (output,err)= p.communicate()
                p_status = p.wait()
 
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--withecalcontainment", action="store_true", help="include ecal containment")
    parser.add_argument("--varyfdv", action="store_true", help=" are we varying the fiducial volume and TPC intsrumented volume")
    parser.add_argument("--varyenergy", action="store_true", help=" do we want to vary the energy according to the dune flux between 1 to 5 GeV")
#    parser.add_argument("--varyenergy", action="store_true", help=" do we want to vary the energy according to the dune flux between 1 to 5 GeV")
    parser.add_argument("--varyBfield", action="store_true", help=" do we want to vary the B field")
    parser.add_argument("--tpcfidrad", type=float, action="store", default=160.0)
    parser.add_argument("--tpcinstrumentedrad", type=float, action="store", default=249.45)
    parser.add_argument("--bfield", type=float, action="store", default=0.5)
    args = parser.parse_args() 
    print(args.tpcfidrad)
    ecalinnerrad = args.tpcinstrumentedrad + 28.57
    ecalouterrad = ecalinnerrad + 44.232
    main()
