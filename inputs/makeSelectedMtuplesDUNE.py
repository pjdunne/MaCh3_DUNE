import os, sys
from ROOT import *
import argparse


if __name__ == "__main__":
    
    # Get arguments (ie. SK file name)
    parser = argparse.ArgumentParser(description='Splits an SK file into five files, containing the numu-,nue-, nue1pi- numucc1pi- selected events')
    parser.add_argument('filename', nargs=1)
    args = parser.parse_args()
    skfile = args.filename[0]

    print "Working on file ", skfile

    # Get SK tree
    f_sk = TFile(skfile, "OPEN")
    if not f_sk:
        print "Error: could not open file ", skfile
        sys.exit()
    sk_tree = f_sk.Get("caf")
    if not sk_tree:
        sk_tree = f_sk.Get("h1")

    # Get entries in SK tree
    n_total = sk_tree.GetEntries()
    print n_total, " events in tree."

    # Get the SK filename (without '.root')
    # skfile = skfile[:-5]

    # Make new filenames for files containing only the numu/nue trees
    numu_numu_numuselec_file  = skfile[:-13]+"_numu_x_numu_numuselec.root"
    numu_numu_nueselec_file  = skfile[:-13]+"_numu_x_numu_nueselec.root"
    numu_nue_numuselec_file   = skfile[:-13]+"_numu_x_nue_numuselec.root"
    numu_nue_nueselec_file   = skfile[:-13]+"_numu_x_nue_nueselec.root"
    numu_nutau_numuselec_file = skfile[:-13]+"_numu_x_nutau_numuselec.root"
    numu_nutau_nueselec_file = skfile[:-13]+"_numu_x_nutau_nueselec.root"
    nue_nue_numuselec_file    = skfile[:-13]+"_nue_x_nue_numuselec.root"
    nue_nue_nueselec_file    = skfile[:-13]+"_nue_x_nue_nueselec.root"
    nue_numu_numuselec_file   = skfile[:-13]+"_nue_x_numu_numuselec.root"
    nue_numu_nueselec_file   = skfile[:-13]+"_nue_x_numu_nueselec.root"
    nue_nutau_numuselec_file  = skfile[:-13]+"_nue_x_nutau_numuselec.root"
    nue_nutau_nueselec_file  = skfile[:-13]+"_nue_x_nutau_nueselec.root"
    numubar_numubar_numuselec_file  = skfile[:-13]+"_numubar_x_numubar_numuselec.root"
    numubar_numubar_nueselec_file  = skfile[:-13]+"_numubar_x_numubar_nueselec.root"
    numubar_nuebar_numuselec_file   = skfile[:-13]+"_numubar_x_nuebar_numuselec.root"
    numubar_nuebar_nueselec_file   = skfile[:-13]+"_numubar_x_nuebar_nueselec.root"
    numubar_nutaubar_numuselec_file = skfile[:-13]+"_numubar_x_nutaubar_numuselec.root"
    numubar_nutaubar_nueselec_file = skfile[:-13]+"_numubar_x_nutaubar_nueselec.root"
    nuebar_nuebar_numuselec_file    = skfile[:-13]+"_nuebar_x_nuebar_numuselec.root"
    nuebar_nuebar_nueselec_file    = skfile[:-13]+"_nuebar_x_nuebar_nueselec.root"
    nuebar_numubar_numuselec_file   = skfile[:-13]+"_nuebar_x_numubar_numuselec.root"
    nuebar_numubar_nueselec_file   = skfile[:-13]+"_nuebar_x_numubar_nueselec.root"
    nuebar_nutaubar_numuselec_file  = skfile[:-13]+"_nuebar_x_nutaubar_numuselec.root"
    nuebar_nutaubar_nueselec_file  = skfile[:-13]+"_nuebar_x_nutaubar_nueselec.root"
    print "Making new files: "
    print numu_numu_numuselec_file
    print numu_numu_nueselec_file
    print numu_nue_numuselec_file
    print numu_nue_nueselec_file
    print numu_nutau_numuselec_file
    print numu_nutau_nueselec_file
    print nue_nue_numuselec_file
    print nue_nue_nueselec_file
    print nue_numu_numuselec_file
    print nue_numu_nueselec_file
    print nue_nutau_numuselec_file
    print nue_nutau_nueselec_file
    print numubar_numubar_numuselec_file
    print numubar_numubar_nueselec_file
    print numubar_nuebar_numuselec_file
    print numubar_nuebar_nueselec_file
    print numubar_nutaubar_numuselec_file
    print numubar_nutaubar_nueselec_file
    print nuebar_nuebar_numuselec_file
    print nuebar_nuebar_nueselec_file
    print nuebar_numubar_numuselec_file
    print nuebar_numubar_nueselec_file
    print nuebar_nutaubar_numuselec_file
    print nuebar_nutaubar_nueselec_file
    print "Check: number of events in original SK tree = ", n_total

## true numu_numu events 

    # Make numu_numu_numuselec tree and save to file
    f_numu_numu_numuselec = TFile(numu_numu_numuselec_file,"RECREATE")
    numu_numu_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==14 and sk_tree.cvnnumu > 0.5: #numu_numu_numuselec selection
            numu_numu_numuselec_tree.Fill()

    n_numu_numu_numuselec = numu_numu_numuselec_tree.GetEntries()
    print "Check: number of events in numu_numu_numuselec tree = ", n_numu_numu_numuselec	
    numu_numu_numuselec_tree.Write("caf")
    f_numu_numu_numuselec.Close()
    if n_numu_numu_numuselec == 0:
       os.remove(numu_numu_numuselec_file)

    # Make numu_numu_nueselec tree and save to file
    f_numu_numu_nueselec = TFile(numu_numu_nueselec_file,"RECREATE")
    numu_numu_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==14 and sk_tree.cvnnue > 0.85: #numu_numu_nueselec selection
            numu_numu_nueselec_tree.Fill()

    n_numu_numu_nueselec = numu_numu_nueselec_tree.GetEntries()
    print "Check: number of events in numu_numu_nueselec tree = ", n_numu_numu_nueselec	
    numu_numu_nueselec_tree.Write("caf")
    f_numu_numu_nueselec.Close()
    if n_numu_numu_nueselec == 0:
       os.remove(numu_numu_nueselec_file)

## true numu_nue events


    # Make numu_nue_numuselec tree and save to file
    f_numu_nue_numuselec = TFile(numu_nue_numuselec_file,"RECREATE")
    numu_nue_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==12 and sk_tree.cvnnumu > 0.5: #numu_nue_numuselec selection
            numu_nue_numuselec_tree.Fill()

    n_numu_nue_numuselec = numu_nue_numuselec_tree.GetEntries()
    print "Check: number of events in numu_nue_numuselec tree = ", n_numu_nue_numuselec

    numu_nue_numuselec_tree.Write("caf")
    f_numu_nue_numuselec.Close()
    if n_numu_nue_numuselec == 0:
       os.remove(numu_nue_numuselec_file)


    # Make numu_nue_nueselec tree and save to file
    f_numu_nue_nueselec = TFile(numu_nue_nueselec_file,"RECREATE")
    numu_nue_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==12 and sk_tree.cvnnue > 0.85: #numu_nue_nueselec selection
            numu_nue_nueselec_tree.Fill()

    n_numu_nue_nueselec = numu_nue_nueselec_tree.GetEntries()
    print "Check: number of events in numu_nue_nueselec tree = ", n_numu_nue_nueselec

    numu_nue_nueselec_tree.Write("caf")
    f_numu_nue_nueselec.Close()
    if n_numu_nue_nueselec == 0:
       os.remove(numu_nue_nueselec_file)

    
## true numu_nutau events


    # Make numu_nutau_numuselec tree and save to file
    f_numu_nutau_numuselec = TFile(numu_nutau_numuselec_file,"RECREATE")
    numu_nutau_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==16 and sk_tree.cvnnumu > 0.5: #numu_nutau_numuselec selection
            numu_nutau_numuselec_tree.Fill()

    n_numu_nutau_numuselec = numu_nutau_numuselec_tree.GetEntries()
    print "Check: number of events in numu_nutau_numuselec tree = ", n_numu_nutau_numuselec

    numu_nutau_numuselec_tree.Write("caf")
    f_numu_nutau_numuselec.Close()
    if n_numu_nutau_numuselec == 0:
       os.remove(numu_nutau_numuselec_file)
    
    

    # Make numu_nutau_nueselec tree and save to file
    f_numu_nutau_nueselec = TFile(numu_nutau_nueselec_file,"RECREATE")
    numu_nutau_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==16 and sk_tree.cvnnue > 0.85: #numu_nutau_nueselec selection
            numu_nutau_nueselec_tree.Fill()

    n_numu_nutau_nueselec = numu_nutau_nueselec_tree.GetEntries()
    print "Check: number of events in numu_nutau_nueselec tree = ", n_numu_nutau_nueselec

    numu_nutau_nueselec_tree.Write("caf")
    f_numu_nutau_nueselec.Close()
    if n_numu_nutau_nueselec == 0:
       os.remove(numu_nutau_nueselec_file)

## true nue_nue events

    
    # Make nue_nue_numuselec tree and save to file
    f_nue_nue_numuselec = TFile(nue_nue_numuselec_file,"RECREATE")
    nue_nue_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==12 and sk_tree.cvnnumu > 0.5: #nue_nue_numuselec selection
            nue_nue_numuselec_tree.Fill()

    n_nue_nue_numuselec = nue_nue_numuselec_tree.GetEntries()
    print "Check: number of events in nue_nue_numuselec tree = ", n_nue_nue_numuselec

    nue_nue_numuselec_tree.Write("caf")
    f_nue_nue_numuselec.Close()
    if n_nue_nue_numuselec == 0:
       os.remove(nue_nue_numuselec_file)

    
    # Make nue_nue_nueselec tree and save to file
    f_nue_nue_nueselec = TFile(nue_nue_nueselec_file,"RECREATE")
    nue_nue_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==12 and sk_tree.cvnnue > 0.85: #nue_nue_nueselec selection
            nue_nue_nueselec_tree.Fill()

    n_nue_nue_nueselec = nue_nue_nueselec_tree.GetEntries()
    print "Check: number of events in nue_nue_nueselec tree = ", n_nue_nue_nueselec

    nue_nue_nueselec_tree.Write("caf")
    f_nue_nue_nueselec.Close()
    if n_nue_nue_nueselec == 0:
       os.remove(nue_nue_nueselec_file)

    
## true nue_numu events


    # Make nue_numu_numuselec tree and save to file
    f_nue_numu_numuselec = TFile(nue_numu_numuselec_file,"RECREATE")
    nue_numu_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==14 and sk_tree.cvnnumu > 0.5: #nue_numu_numuselec selection
            nue_numu_numuselec_tree.Fill()

    n_nue_numu_numuselec = nue_numu_numuselec_tree.GetEntries()
    print "Check: number of events in nue_numu_numuselec tree = ", n_nue_numu_numuselec

    nue_numu_numuselec_tree.Write("caf")
    f_nue_numu_numuselec.Close()
    if n_nue_numu_numuselec == 0:
       os.remove(nue_numu_numuselec_file)

    # Make nue_numu_nueselec tree and save to file
    f_nue_numu_nueselec = TFile(nue_numu_nueselec_file,"RECREATE")
    nue_numu_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==14 and sk_tree.cvnnue > 0.85: #nue_numu_nueselec selection
            nue_numu_nueselec_tree.Fill()

    n_nue_numu_nueselec = nue_numu_nueselec_tree.GetEntries()
    print "Check: number of events in nue_numu_nueselec tree = ", n_nue_numu_nueselec

    nue_numu_nueselec_tree.Write("caf")
    f_nue_numu_nueselec.Close()
    if n_nue_numu_nueselec == 0:
       os.remove(nue_numu_nueselec_file)
    

## true nue_nutau events


    # Make nue_nutau_numuselec tree and save to file
    f_nue_nutau_numuselec = TFile(nue_nutau_numuselec_file,"RECREATE")
    nue_nutau_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==16 and sk_tree.cvnnumu > 0.5: #nue_nutau_numuselec
            nue_nutau_numuselec_tree.Fill()

    n_nue_nutau_numuselec = nue_nutau_numuselec_tree.GetEntries()
    print "Check: number of events in nue_nutau_numuselec tree = ", n_nue_nutau_numuselec

    nue_nutau_numuselec_tree.Write("caf")
    f_nue_nutau_numuselec.Close()
    if n_nue_nutau_numuselec == 0:
       os.remove(nue_nutau_numuselec_file)
    
    
    # Make nue_nutau_nueselec tree and save to file
    f_nue_nutau_nueselec = TFile(nue_nutau_nueselec_file,"RECREATE")
    nue_nutau_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==16 and sk_tree.cvnnue > 0.85: #nue_nutau_nueselec
            nue_nutau_nueselec_tree.Fill()

    n_nue_nutau_nueselec = nue_nutau_nueselec_tree.GetEntries()
    print "Check: number of events in nue_nutau_nueselec tree = ", n_nue_nutau_nueselec

    nue_nutau_nueselec_tree.Write("caf")
    f_nue_nutau_nueselec.Close()
    if n_nue_nutau_nueselec == 0:
       os.remove(nue_nutau_nueselec_file)
    
## true numubar_numubar events


    # Make numubar_numubar_numuselec tree and save to file
    f_numubar_numubar_numuselec = TFile(numubar_numubar_numuselec_file,"RECREATE")
    numubar_numubar_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-14 and sk_tree.cvnnumu > 0.5: #numubar_numubar_numuselec
            numubar_numubar_numuselec_tree.Fill()

    n_numubar_numubar_numuselec = numubar_numubar_numuselec_tree.GetEntries()
    print "Check: number of events in numubar_numubar_numuselec tree = ", n_numubar_numubar_numuselec

    numubar_numubar_numuselec_tree.Write("caf")
    f_numubar_numubar_numuselec.Close()
    if n_numubar_numubar_numuselec == 0:
       os.remove(numubar_numubar_numuselec_file)

    # Make numubar_numubar_nueselec tree and save to file
    f_numubar_numubar_nueselec = TFile(numubar_numubar_nueselec_file,"RECREATE")
    numubar_numubar_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-14 and sk_tree.cvnnue > 0.85: #numubar_numubar_nueselec
            numubar_numubar_nueselec_tree.Fill()

    n_numubar_numubar_nueselec = numubar_numubar_nueselec_tree.GetEntries()
    print "Check: number of events in numubar_numubar_nueselec tree = ", n_numubar_numubar_nueselec

    numubar_numubar_nueselec_tree.Write("caf")
    f_numubar_numubar_nueselec.Close()
    if n_numubar_numubar_nueselec == 0:
       os.remove(numubar_numubar_nueselec_file)
    
    
## true numubar_nuebar events


    # Make numubar_nuebar_numuselec tree and save to file
    f_numubar_nuebar_numuselec = TFile(numubar_nuebar_numuselec_file,"RECREATE")
    numubar_nuebar_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-12 and sk_tree.cvnnumu > 0.5: #numubar_nuebar_numuselec
            numubar_nuebar_numuselec_tree.Fill()

    n_numubar_nuebar_numuselec = numubar_nuebar_numuselec_tree.GetEntries()
    print "Check: number of events in numubar_nuebar_numuselec tree = ", n_numubar_nuebar_numuselec

    numubar_nuebar_numuselec_tree.Write("caf")
    f_numubar_nuebar_numuselec.Close()
    if n_numubar_nuebar_numuselec == 0:
       os.remove(numubar_nuebar_numuselec_file)

    # Make numubar_nuebar_nueselec tree and save to file
    f_numubar_nuebar_nueselec = TFile(numubar_nuebar_nueselec_file,"RECREATE")
    numubar_nuebar_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-12 and sk_tree.cvnnue > 0.85: #numubar_nuebar_nueselec
            numubar_nuebar_nueselec_tree.Fill()

    n_numubar_nuebar_nueselec = numubar_nuebar_nueselec_tree.GetEntries()
    print "Check: number of events in numubar_nuebar_nueselec tree = ", n_numubar_nuebar_nueselec

    numubar_nuebar_nueselec_tree.Write("caf")
    f_numubar_nuebar_nueselec.Close()
    if n_numubar_nuebar_nueselec == 0:
       os.remove(numubar_nuebar_nueselec_file)
    


## true numubar_nutaubar events


    # Make numubar_nutaubar_numuselec tree and save to file
    f_numubar_nutaubar_numuselec = TFile(numubar_nutaubar_numuselec_file,"RECREATE")
    numubar_nutaubar_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-16 and sk_tree.cvnnumu > 0.5: #numubar_nutaubar_numuselec
            numubar_nutaubar_numuselec_tree.Fill()

    n_numubar_nutaubar_numuselec = numubar_nutaubar_numuselec_tree.GetEntries()
    print "Check: number of events in numubar_nutaubar_numuselec tree = ", n_numubar_nutaubar_numuselec

    numubar_nutaubar_numuselec_tree.Write("caf")
    f_numubar_nutaubar_numuselec.Close()
    if n_numubar_nutaubar_numuselec == 0:
       os.remove(numubar_nutaubar_numuselec_file)

    
    # Make numubar_nutaubar_nueselec tree and save to file
    f_numubar_nutaubar_nueselec = TFile(numubar_nutaubar_nueselec_file,"RECREATE")
    numubar_nutaubar_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-16 and sk_tree.cvnnue > 0.85: #numubar_nutaubar_nueselec
            numubar_nutaubar_nueselec_tree.Fill()

    n_numubar_nutaubar_nueselec = numubar_nutaubar_nueselec_tree.GetEntries()
    print "Check: number of events in numubar_nutaubar_nueselec tree = ", n_numubar_nutaubar_nueselec

    numubar_nutaubar_nueselec_tree.Write("caf")
    f_numubar_nutaubar_nueselec.Close()
    if n_numubar_nutaubar_nueselec == 0:
       os.remove(numubar_nutaubar_nueselec_file)

## true nuebar_nuebar events
    
    
    # Make nuebar_nuebar_numuselec tree and save to file
    f_nuebar_nuebar_numuselec = TFile(nuebar_nuebar_numuselec_file,"RECREATE")
    nuebar_nuebar_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-12 and sk_tree.cvnnumu > 0.5: #nuebar_nuebar_numuselec
            nuebar_nuebar_numuselec_tree.Fill()

    n_nuebar_nuebar_numuselec = nuebar_nuebar_numuselec_tree.GetEntries()
    print "Check: number of events in nuebar_nuebar_numuselec tree = ", n_nuebar_nuebar_numuselec

    nuebar_nuebar_numuselec_tree.Write("caf")
    f_nuebar_nuebar_numuselec.Close()
    if n_nuebar_nuebar_numuselec == 0:
       os.remove(nuebar_nuebar_numuselec_file)

    # Make nuebar_nuebar_nueselec tree and save to file
    f_nuebar_nuebar_nueselec = TFile(nuebar_nuebar_nueselec_file,"RECREATE")
    nuebar_nuebar_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-12 and sk_tree.cvnnue > 0.85: #nuebar_nuebar_nueselec
            nuebar_nuebar_nueselec_tree.Fill()

    n_nuebar_nuebar_nueselec = nuebar_nuebar_nueselec_tree.GetEntries()
    print "Check: number of events in nuebar_nuebar_nueselec tree = ", n_nuebar_nuebar_nueselec

    nuebar_nuebar_nueselec_tree.Write("caf")
    f_nuebar_nuebar_nueselec.Close()
    if n_nuebar_nuebar_nueselec == 0:
       os.remove(nuebar_nuebar_nueselec_file)
    

## true nuebar_numubar events


    # Make nuebar_numubar_numuselec tree and save to file
    f_nuebar_numubar_numuselec = TFile(nuebar_numubar_numuselec_file,"RECREATE")
    nuebar_numubar_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-14 and sk_tree.cvnnumu > 0.5: #nuebar_numubar_numuselec
            nuebar_numubar_numuselec_tree.Fill()

    n_nuebar_numubar_numuselec = nuebar_numubar_numuselec_tree.GetEntries()
    print "Check: number of events in nuebar_numubar_numuselec tree = ", n_nuebar_numubar_numuselec

    nuebar_numubar_numuselec_tree.Write("caf")
    f_nuebar_numubar_numuselec.Close()
    if n_nuebar_numubar_numuselec == 0:
       os.remove(nuebar_numubar_numuselec_file)

    # Make nuebar_numubar_nueselec tree and save to file
    f_nuebar_numubar_nueselec = TFile(nuebar_numubar_nueselec_file,"RECREATE")
    nuebar_numubar_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-14 and sk_tree.cvnnue > 0.85: #nuebar_numubar_nueselec
            nuebar_numubar_nueselec_tree.Fill()

    n_nuebar_numubar_nueselec = nuebar_numubar_nueselec_tree.GetEntries()
    print "Check: number of events in nuebar_numubar_nueselec tree = ", n_nuebar_numubar_nueselec

    nuebar_numubar_nueselec_tree.Write("caf")
    f_nuebar_numubar_nueselec.Close()
    if n_nuebar_numubar_nueselec == 0:
       os.remove(nuebar_numubar_nueselec_file)
    

## true nuebar_nutaubar events


    # Make nuebar_nutaubar_numuselec tree and save to file
    f_nuebar_nutaubar_numuselec = TFile(nuebar_nutaubar_numuselec_file,"RECREATE")
    nuebar_nutaubar_numuselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-16 and sk_tree.cvnnumu > 0.5: #nuebar_nutaubar_numuselec
            nuebar_nutaubar_numuselec_tree.Fill()

    n_nuebar_nutaubar_numuselec = nuebar_nutaubar_numuselec_tree.GetEntries()
    print "Check: number of events in nuebar_nutaubar_numuselec tree = ", n_nuebar_nutaubar_numuselec

    nuebar_nutaubar_numuselec_tree.Write("caf")
    f_nuebar_nutaubar_numuselec.Close()
    if n_nuebar_nutaubar_numuselec == 0:
       os.remove(nuebar_nutaubar_numuselec_file)
    

    # Make nuebar_nutaubar_nueselec tree and save to file
    f_nuebar_nutaubar_nueselec = TFile(nuebar_nutaubar_nueselec_file,"RECREATE")
    nuebar_nutaubar_nueselec_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-16 and sk_tree.cvnnue > 0.85: #nuebar_nutaubar_nueselec
            nuebar_nutaubar_nueselec_tree.Fill()

    n_nuebar_nutaubar_nueselec = nuebar_nutaubar_nueselec_tree.GetEntries()
    print "Check: number of events in nuebar_nutaubar_nueselec tree = ", n_nuebar_nutaubar_nueselec

    nuebar_nutaubar_nueselec_tree.Write("caf")
    f_nuebar_nutaubar_nueselec.Close()
    if n_nuebar_nutaubar_nueselec == 0:
       os.remove(nuebar_nutaubar_nueselec_file)


    print "Sum of numuselec event selections = ", n_numu_numu_numuselec  + n_numu_nue_numuselec + n_numu_nutau_numuselec + n_nue_nue_numuselec + n_nue_numu_numuselec + n_nue_nutau_numuselec + n_numubar_numubar_numuselec + n_numubar_nuebar_numuselec + n_numubar_nutaubar_numuselec +  n_nuebar_nuebar_numuselec + n_nuebar_numubar_numuselec + n_nuebar_nutaubar_numuselec


    print "Sum of nueselec event selections = ", n_numu_numu_nueselec  + n_numu_nue_nueselec + n_numu_nutau_nueselec + n_nue_nue_nueselec + n_nue_numu_nueselec + n_nue_nutau_nueselec + n_numubar_numubar_nueselec + n_numubar_nuebar_nueselec + n_numubar_nutaubar_nueselec +  n_nuebar_nuebar_nueselec + n_nuebar_numubar_nueselec + n_nuebar_nutaubar_nueselec


    print "Sum of all  event selections = ", n_numu_numu_numuselec  + n_numu_nue_numuselec + n_numu_nutau_numuselec + n_nue_nue_numuselec + n_nue_numu_numuselec + n_nue_nutau_numuselec + n_numubar_numubar_numuselec + n_numubar_nuebar_numuselec + n_numubar_nutaubar_numuselec +  n_nuebar_nuebar_numuselec + n_nuebar_numubar_numuselec + n_nuebar_nutaubar_numuselec + n_numu_numu_nueselec  + n_numu_nue_nueselec + n_numu_nutau_nueselec + n_nue_nue_nueselec + n_nue_numu_nueselec + n_nue_nutau_nueselec + n_numubar_numubar_nueselec + n_numubar_nuebar_nueselec + n_numubar_nutaubar_nueselec +  n_nuebar_nuebar_nueselec + n_nuebar_numubar_nueselec + n_nuebar_nutaubar_nueselec



    print "This should be equal to the total number of events: ", n_total 
