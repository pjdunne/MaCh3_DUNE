import subprocess
import os
from ruamel.yaml import YAML
import argparse

def modify_samplepdfconfig(config_file_path, isaccepted=[0.0, 2.0], enubounds=[2.0, 3.0], tpcfidrad=160.0, tpcinstrumentedrad=249.45, ecalinnerrad=278.02, ecalouterrad=322.252, adcsamplingfreq=20.0, spatialresolution=2.5, bfield=0.5, KinematicStr=True):
    yaml=YAML()
    with open(config_file_path) as f:
        conf=yaml.load(f)

# Find the index of the line after the search string
    found = False
    index_to_modify =0
    if KinematicStr:
        for obj in conf['SelectionCuts']:
#       print(obj)
            if obj['KinematicStr'] == 'TrueNeutrinoEnergy':
 #          print(len(obj))
                obj['Bounds'] = enubounds
 #              print(newstring)
            elif obj['KinematicStr'] == 'IsAccepted':
                obj['Bounds'] = isaccepted
                
    if args.withecalcontainment:
        conf['SampleBools']['ecal_containment'] = 'yes'
    elif not args.withecalcontainment:
        conf['SampleBools']['ecal_containment'] = 'no'
    conf['SampleCuts']['TPCFidRadius'] = tpcfidrad
    conf['SampleCuts']['TPCInstrumentedRadius'] = tpcinstrumentedrad
    conf['SampleCuts']['ECALInnerRadius'] = ecalinnerrad
    conf['SampleCuts']['ECALOuterRadius'] = ecalouterrad
    conf['SampleCuts']['B_field']=bfield
    conf['SampleCuts']['adc_sampling_frequency']=adcsamplingfreq
    conf['SampleCuts']['spatial_resolution']=spatialresolution
    with open(config_file_path, 'w') as f:
        yaml.dump(conf, f)

def modify_config_selections(config_file_path, newstring):
    yaml=YAML()
    with open(config_file_path) as f:
        conf=yaml.load(f)
    conf['General']['Output']['FileName'] = newstring
    print(conf)

    with open(config_file_path, 'w') as f:
        yaml.dump(conf, f)
        
def main():
    if args.withecalcontainment:
        ecalcontainment = "withecalcontainment"
    else:
        ecalcontainment = "noecalcontainment"
    if args.varyspatialresolution:
        spatialres = [1.0, 2.0, 2.5, 3.0]
    else:
        spatialres = [args.spatialresolution]
    if args.varyadcfrequency:
        adcfreq = [2.0, 10.0, 20.0]
    else:
        adcfreq = [args.adcsamplingfrequency]
    if args.varyfdv:
        fdv  = [30.0, 50.0, 70.0, 90.0, 110.0]
    else:
        fdv = [args.tpcfidrad]
    if args.varytpcinstrumentedrad:
        instrumentedrad = [199.45, 249.45]
        ecalinnerrad = [228.02, 278.02]
        ecalouterrad = [272.252, 322.252]
    else:
        instrumentedrad = [args.tpcinstrumentedrad] 
    if args.varyenergy:
        enu = [0.00, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 4.00, 4.50, 5.00]
    else:
        enu = [2.0, 3.0]
    if args.varyBfield:
        bfield = [0.5, 0.4, 0.3, 0.2, 0.1]
    else:
        bfield = [0.5]
    print(instrumentedrad)
    for i_instrumentedrad in range(0, len(instrumentedrad)):
        for i_fdv in range(0, len(fdv)):
            fidrad = round(instrumentedrad[i_instrumentedrad]-fdv[i_fdv], 3)
            if(fidrad>=50.0):
                for i_bfield in range(0, len(bfield)):
                    for i_res in range(0, len(spatialres)):
                        for i_freq in range(0, len(adcfreq)):
#                   enu = [0.00, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 4.00, 4.50, 5.00]
#                   enu=[1.50, 1.75]
                            accepted_arr=[0.0, 1.0, 2.0]
                            k=0
                            print(len(accepted_arr))
#                            newdir='configs/SamplePDFConfigs_{pointres}mmpointres_{samplingfreq}hz_fidrad{fidrad}_Bfield{bfield}_instrumentedrad{instrumentedradius}_{ecal_containment}'.format(pointres=str(spatialres[i_res]).replace(".", "_"), samplingfreq=str(adcfreq[i_freq]).replace(".", "_"),fidrad=f"{fidrad:.2f}".replace(".", "_"), bfield = str(bfield[i_bfield]).replace(".", "_"), instrumentedradius = str(instrumentedrad[i_instrumentedrad]).replace(".", "_"), ecal_containment=ecalcontainment)
                            outdir='outputs/{pointres}mmpointres_{samplingfreq}hz_fidrad{fidrad}_Bfield{bfield}_instrumentedrad{instrumentedradius}_{ecal_containment}'.format(pointres=str(spatialres[i_res]).replace(".", "_"), samplingfreq=str(adcfreq[i_freq]).replace(".", "_"), fidrad=f"{fidrad:.2f}".replace(".", "_"), bfield = str(bfield[i_bfield]).replace(".", "_"), instrumentedradius = str(instrumentedrad[i_instrumentedrad]).replace(".", "_"), ecal_containment=ecalcontainment)
                           # p_1 = subprocess.Popen('mkdir {newdir}', shell=True)
#                            p_1 = subprocess.Popen(['mkdir', '-p', f"{newdir}"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#                            (output,err)= p_1.communicate()
#                            p_1status = p_1.wait()       
                            p_out = subprocess.Popen(['mkdir', '-p', f"{outdir}"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            (output,err)= p_out.communicate()
                            p_outstatus = p_out.wait()       

                            for i in range(0, len(enu)-1):
#                                config_samplepdf = 'configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml'
                                #kinematicstr = 'TrueNeutrinoEnergy'
                    #            newvalues = '[{boundlow:.2f}, {boundup:.2f}]'.format(boundlow=enu[i], boundup=enu[i+1])
                                #newvalues = [enu[i], enu[i+1]]
                                #modify_config(config_samplepdf, 'SelectionCuts', kinematicstr, newvalues, True)
                    #            modify_config(config_samplepdf, 'SelectionCuts', kinematicstr, newvalues, True)
                                for j in range(0, len(accepted_arr)-1):
                                    print(accepted_arr)
                                    print(j)
 #                                   config_samplepdf_new = '{newdir}/SamplePDFDuneNDGAr_FHC_CCnumuselec_{k}.yaml'.format(newdir=newdir, k=k)
 #                                   p_2 = subprocess.Popen(['cp', f"{config_samplepdf}", f"{config_samplepdf_new}"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
##                                    (output,err)= p_2.communicate()
#                                    p_2status = p_2.wait()
#                                    modify_samplepdfconfig(config_samplepdf_new, [accepted_arr[j], 2.0], enubounds=[enu[i], enu[i+1]], tpcfidrad=fidrad, tpcinstrumentedrad=instrumentedrad[i_instrumentedrad], ecalinnerrad=ecalinnerrad[i_instrumentedrad], ecalouterrad=ecalouterrad[i_instrumentedrad], adcsamplingfreq=adcfreq[i_freq], spatialresolution=spatialres[i_res], bfield=bfield[i_bfield], KinematicStr=True)
#                                    #modify_config(config_samplepdf_new, 'SelectionCuts', 'IsAccepted', [accepted_arr[j], 2.0], True) 
#                                    config_selec = 'configs/Selections_NDGAr_acceptancecorrection_2dhist.yaml'
#                                    config_selec_new = '{newdir}/Selections_NDGAr_acceptancecorrection_2dhist_{k}.yaml'.format(newdir=newdir,k=k)
 #                                   p_3 = subprocess.Popen(['cp', f"{config_selec}", f"{config_selec_new}"],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
 #                                   (output,err)= p_3.communicate()
 #                                   p_3status = p_3.wait()
 #                                   nameout = 'Output: '
 #                                   if j == 1:
 #                                       accepted = "Accepted"
  #                                  else:
  #                                      accepted = "All"
                                   # filename = 'Selections_CC_AcceptanceCorrection_{accepted}Events_Enubin{binnum}'.format(accepted=accepted, binnum=i)
#                                    filename = '{outdir}/Selections_CC_AcceptanceCorrection_{accepted}Events_Enubin{binnum}'.format(outdir=outdir, accepted=accepted, binnum=i)
#                                    newoutput = '{filename}.root'.format(filename=filename)
#                                    modify_config_selections(config_selec_new, newoutput)
                    #                file = 'acceptance_{inum}_{jnum}.job'.format(inum=i, jnum=j)
                    #                command = 'Selections {config_selec_new} {config_samplepdf_new} &> {filename}.txt'.format(filename=filename)
                    #                print(command)
                    #                p = subprocess.Popen(command, shell=True)
                    #                (output,err)= p.communicate()
                    #                p_status = p.wait()
                                    k=k+1
 
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--withecalcontainment", action="store_true", help="include ecal containment")
    parser.add_argument("--varyfdv", action="store_true", help=" are we varying the fiducial volume")
    parser.add_argument("--varytpcinstrumentedrad", action="store_true", help=" are we varying the tpc instrumented radius")
    parser.add_argument("--varyenergy", action="store_true", help=" do we want to vary the energy according to the dune flux between 1 to 5 GeV")
    parser.add_argument("--varyspatialresolution", action="store_true", help=" do we want to vary the spatial resolution in the readout plane")
    parser.add_argument("--varyadcfrequency", action="store_true", help=" do we want to vary the adc sampling frequency")
    parser.add_argument("--varyBfield", action="store_true", help=" do we want to vary the B field")
    parser.add_argument("--tpcfidrad", type=float, action="store", default=160.0)
    parser.add_argument("--tpcinstrumentedrad", type=float, action="store", default=249.45)
    parser.add_argument("--bfield", type=float, action="store", default=0.5)
    parser.add_argument("--spatialresolution", type=float, action="store", default=2.5)
    parser.add_argument("--adcsamplingfrequency", type=float, action="store", default=60.0)

    args = parser.parse_args() 
    print(args.tpcfidrad)
    ecalinnerrad = round(args.tpcinstrumentedrad + 28.57, 3)
    ecalouterrad = round(ecalinnerrad + 44.232, 3)
    main()
