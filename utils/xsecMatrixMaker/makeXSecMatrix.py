#!/usr/bin/python

import xml.etree.ElementTree as ET
import ROOT
import math
import sys

# LOL DOCUMENTATION
# ACK, Oct 2016
#
# Produces a MaCh3 compatible matrix for cross section systematics
# First argument is the xml file
# Second argument is the output root file
# Code is hugely hacked and ugly, but totes functional
#

if len(sys.argv) != 3:
  print "Sorry, I need two arguments"
  print "./makeXSecMatrix.py input.xml output.root"
  sys.exit()

tree = ET.parse(sys.argv[1])
fromxml = tree.getroot()

maxelements=300
nummodes=12


xsec_param_names        = ROOT.TObjArray()
xsec_param_prior_unnorm = ROOT.TVectorD(maxelements)
xsec_param_nom_unnorm   = ROOT.TVectorD(maxelements)
xsec_param_prior        = ROOT.TVectorD(maxelements)
xsec_param_nom          = ROOT.TVectorD(maxelements)
xsec_param_lb           = ROOT.TVectorD(maxelements)
xsec_param_ub           = ROOT.TVectorD(maxelements)
xsec_error              = ROOT.TVectorD(maxelements)
xsec_stepscale          = ROOT.TVectorD(maxelements)
xsec_param_id           = ROOT.TMatrixD(maxelements, 2)

xsec_sk_spline_names    = ROOT.TObjArray()
xsec_sk_spline_modes    = ROOT.TObjArray()

xsec_nd_spline_names    = ROOT.TObjArray()

xsec_norm_etru_bnd_low  = ROOT.TObjArray()
xsec_norm_etru_bnd_high = ROOT.TObjArray()
xsec_norm_q2_true_bnd_low  = ROOT.TObjArray()
xsec_norm_q2_true_bnd_high = ROOT.TObjArray()
xsec_norm_modes         = ROOT.TObjArray()
xsec_norm_elements      = ROOT.TObjArray()
xsec_norm_nupdg         = ROOT.TObjArray()
xsec_norm_prod_nupdg         = ROOT.TObjArray()
xsec_norm_horncurrents         = ROOT.TObjArray()

corr_dicts = []
corr_names = []


nelem = 0
for child in fromxml:

    # Get the name attribute
    name = ROOT.TObjString(child.attrib['name'])


    # Add the name to the TObjArray
    xsec_param_names.AddLast(name)
    corr_names.append(name)

    # Get the prior from attribute
    prior = float(child.attrib['prior'])
    xsec_param_prior_unnorm[nelem]=float(prior)

    # Get the nominal from attribute
    nom = float(child.attrib['nom'])
    xsec_param_nom_unnorm[nelem]=float(nom)

    # Get the lower bound from attribute
    lb = float(child.attrib['lb'])
    ub = float(child.attrib['ub'])
    error = float(child.attrib['error'])
    renorm = child.attrib['renorm']

    stepscale = float(child.attrib['stepscale'])
    xsec_stepscale[nelem] = float(stepscale)

    # Read in what to of detector the parameter applies to
    # This is meant to be used wtih binary operator in MaCh3, e.g. 25 means apply both to ND280 and SK, 1 means apply only to ND280, etc
    # Is used because for some analyses we used parameters only for numu samples at SK, which had detid 8
    # Could be switched to simple int comparisons...
    xsec_param_id[nelem][1]=int(child.attrib['detid'])

    # If we want to renormalize so that the NEUT default = 1
    if (int(renorm) == 1):
        xsec_param_nom[nelem] = nom/nom
        xsec_param_prior[nelem] = prior/nom
        xsec_error[nelem] = error/nom
        if (lb != -9999):
            xsec_param_lb[nelem] = lb/nom
        else:
            xsec_param_lb[nelem] = lb
        if (ub!= 9999):
            xsec_param_ub[nelem] = ub/nom
        else:
            xsec_param_ub[nelem] = ub
    else:
        xsec_param_nom[nelem] = nom
        xsec_param_prior[nelem] = prior
        xsec_param_lb[nelem] = lb
        xsec_param_ub[nelem] = ub
        xsec_error[nelem] = error

    # Find the correlations
    correlations = child.findall('correlation')
    dict = {}
    for cor in correlations:
        dict[cor.attrib['par']]=float(cor.text)
    corr_dicts.append(dict)

    # Find what type of parameter we're dealing with
    type = str(child.attrib['type'])


    # Spline parameters (e.g. MAQE)

    #All variables need the spline par entry


    if(type=="spline"):
        xsec_param_id[nelem][0]=int(child.attrib['splineind'])
        # Get the SK spline file name attribute
        print xsec_param_id[nelem][1]

        if 'sk_spline_name' in child.attrib:
          sk_spline_name = ROOT.TObjString(child.attrib['sk_spline_name'])
        else: 
          sk_spline_name=ROOT.TObjString("")
        xsec_sk_spline_names.AddAtAndExpand(sk_spline_name, nelem)
        if "_" in str(sk_spline_name):
          print "ERROR : ", sk_spline_name, " contains _ (this is a problem for MaCh3!)"
          sys.exit()


        # Find the mode which this parameter might apply to
        spline_modes = child.findall('sk_mode')
        if(len(spline_modes)==0):
          print "Found no mode list for ",name," filling with blank which means apply to all"
          dummyvec = ROOT.TVectorD(0)
          xsec_sk_spline_modes.AddAtAndExpand(dummyvec,nelem)
        else:
          for mode in spline_modes:
            if mode.text == None or mode.tail == None:
              print "Found empty mode list for",name," filling with blank which means apply to all"
              mode.text= ""
            # Split multiple entries up
            modelist = mode.text.split()
            # Make a TVectorD containing each mode
            modevec = ROOT.TVectorD(len(modelist))
            for i in range(len(modelist)):
              modevec[i]=int(modelist[i])
            xsec_sk_spline_modes.AddAtAndExpand(modevec, nelem)

        if 'nd_spline_name' in child.attrib:
          nd_spline_name = ROOT.TObjString(child.attrib['nd_spline_name'])
        else: 
          nd_spline_name=ROOT.TObjString("")
        xsec_nd_spline_names.AddAtAndExpand(nd_spline_name, nelem)


    #End Spline pars

    # Normalisation parameter (e.g. 2p2h neutrino normalisation)
    elif(type=="norm"):
        xsec_param_id[nelem][0]=-1
    # Functional parameter (e.g. BeRPA)
    elif(type=="function"):
        xsec_param_id[nelem][0]=-2
    else:
        print "Wrong parameter type!"
        quit()

    # Find the mode which this parameter might apply to
    normmodes = child.findall('mode')
    if(len(normmodes)==0 and type=="norm"):
        print "Found no mode list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_modes.AddAtAndExpand(dummyvec,nelem)
    else:
      for mode in normmodes:
        if mode.text == None or mode.tail == None:
          print "Found empty mode list for",name," filling with blank which means apply to all"
          mode.text= ""
        # Split multiple entries up
        modelist = mode.text.split()
        # Make a TVectorD containing each mode
        modevec = ROOT.TVectorD(len(modelist))
        for i in range(len(modelist)):
            modevec[i]=int(modelist[i])
        xsec_norm_modes.AddAtAndExpand(modevec, nelem)

    # Find the target element
    normelements = child.findall('element')
    if(len(normelements)==0 and type=="norm"):
        print "Found no target list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_elements.AddAtAndExpand(dummyvec,nelem)
    else:
      for elem in normelements:
        if  elem.text == None or elem.tail == None:
          print "Found empty target list for ",name," filling with blank which means apply to all"
          elem.text = ""
        elemlist = elem.text.split()
        elemvec = ROOT.TVectorD(len(elemlist))
        for i in range(len(elemlist)):
            elemvec[i]=int(elemlist[i])
        xsec_norm_elements.AddAtAndExpand(elemvec, nelem)

    # Find the mode which this parameter might apply to
    horncurrents = child.findall('horncurrent')
    if(len(horncurrents)==0 and type=="norm"):
        print "Found no horn current list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_horncurrents.AddAtAndExpand(dummyvec,nelem)
    else:
      for horncurrent in horncurrents:
        if horncurrent.text == None or horncurrent.tail == None:
          print "Found empty horncurrent list for",name," filling with blank which means apply to all"
          horncurrent.text= ""
        # Split multiple entries up
        horncurrentlist = horncurrent.text.split()
        # Make a TVectorD containing each mode
        horncurrentvec = ROOT.TVectorD(len(horncurrentlist))
        for i in range(len(horncurrentlist)):
            horncurrentvec[i]=int(horncurrentlist[i])
        xsec_norm_horncurrents.AddAtAndExpand(horncurrentvec, nelem)

    # Find the neutrino type
    normpdg = child.findall('nupdg')
    if(len(normpdg)==0 and type=="norm"):
        print "Found no nupdg list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_nupdg.AddAtAndExpand(dummyvec,nelem)
    else:
      for pdg in normpdg:
        if pdg.text == None or pdg.tail == None:
          print "Found empty nupdg list for ",name," filling with blank which means apply to all"
          pdg.text = ""
        pdglist = pdg.text.split()
        pdgvec = ROOT.TVectorD(len(pdglist))
        for i in range(len(pdglist)):
            pdgvec[i]=int(pdglist[i])
        xsec_norm_nupdg.AddAtAndExpand(pdgvec, nelem)

    normprodpdg = child.findall('prod_nupdg')
    if(len(normprodpdg)==0 and type=="norm"):
        print "Found no prod_nupdg list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_prod_nupdg.AddAtAndExpand(dummyvec,nelem)
    else:
      for prodpdg in normprodpdg:
        if prodpdg.text == None or prodpdg.tail == None:
          print "Found empty prod_nupdg list for ",name," filling with blank which means apply to all"
          prodpdg.text = ""
        prodpdglist = prodpdg.text.split()
        prodpdgvec = ROOT.TVectorD(len(prodpdglist))
        for i in range(len(prodpdglist)):
            prodpdgvec[i]=int(prodpdglist[i])
        xsec_norm_prod_nupdg.AddAtAndExpand(prodpdgvec, nelem)

    normetrubndlow = child.findall('etru_bnd_low')
    if(len(normetrubndlow)==0 and type=="norm"):
        print "Found no etrubndslow list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(1)
        dummyvec[0]=-999
        xsec_norm_etru_bnd_low.AddAtAndExpand(dummyvec,nelem)
    else:
      for etrubndslow in normetrubndlow:
        if etrubndslow.text == None or etrubndslow.tail == None:
          print "Found empty etrubndslow list for ",name," filling with -999 which means don't apply a bound"
          etrubndslow.text = "-999"
        # Split multiple entries up
        etrubndlowlist = etrubndslow.text.split()
        # Make a TVectorD containing each mode
        etrubndlowvec = ROOT.TVectorD(len(etrubndlowlist))
        for i in range(len(etrubndlowlist)):
            etrubndlowvec[i]=float(etrubndlowlist[i])
        xsec_norm_etru_bnd_low.AddAtAndExpand(etrubndlowvec, nelem)

    normetrubndhigh = child.findall('etru_bnd_high')
    if(len(normetrubndhigh)==0 and type=="norm"):
        print "Found no etrubndshigh list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(1)
        dummyvec[0]=-999
        xsec_norm_etru_bnd_high.AddAtAndExpand(dummyvec,nelem)
    else:
      for etrubndshigh in normetrubndhigh:
        if etrubndshigh.text == None or etrubndshigh.tail == None:
          print "Found empty etrubndshigh list for ",name," filling with -999 which means don't apply a bound"
          etrubndshigh.text = "-999"
        # Split multiple entries up
        etrubndhighlist = etrubndshigh.text.split()
        # Make a TVectorD containing each mode
        etrubndhighvec = ROOT.TVectorD(len(etrubndhighlist))
        for i in range(len(etrubndhighlist)):
            etrubndhighvec[i]=float(etrubndhighlist[i])
        xsec_norm_etru_bnd_high.AddAtAndExpand(etrubndhighvec, nelem)


    normq2truebndlow = child.findall('q2_true_bnd_low')
    if(len(normq2truebndlow)==0 and type=="norm"):
        print "Found no q2_true_bnd_low list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(1)
        dummyvec[0]=-999
        xsec_norm_q2_true_bnd_low.AddAtAndExpand(dummyvec,nelem)
    else:
      for q2_truebndslow in normq2truebndlow:
        if q2_truebndslow.text == None or q2_truebndslow.tail == None:
          print "Found empty q2_truebndslow list for ",name," filling with -999 which means don't apply a bound"
          q2_truebndslow.text = "-999"
        # Split multiple entries up
        q2truebndlowlist = q2_truebndslow.text.split()
        # Make a TVectorD containing each mode
        q2truebndlowvec = ROOT.TVectorD(len(q2truebndlowlist))
        for i in range(len(q2truebndlowlist)):
            q2truebndlowvec[i]=float(q2truebndlowlist[i])
        xsec_norm_q2_true_bnd_low.AddAtAndExpand(q2truebndlowvec, nelem)

    normq2truebndhigh = child.findall('q2_true_bnd_high')
    if(len(normq2truebndhigh)==0 and type=="norm"):
        print "Found no q2_true_bnd_high list for ",name," filling with blank which means apply to all"
        dummyvec = ROOT.TVectorD(1)
        dummyvec[0]=-999
        xsec_norm_q2_true_bnd_high.AddAtAndExpand(dummyvec,nelem)
    else:
      for q2_truebndshigh in normq2truebndhigh:
        if q2_truebndshigh.text == None or q2_truebndshigh.tail == None:
          print "Found empty q2_truebndshigh list for ",name," filling with -999 which means don't apply a bound"
          q2_truebndshigh.text = "-999"
        # Split multiple entries up
        q2truebndhighlist = q2_truebndshigh.text.split()
        # Make a TVectorD containing each mode
        q2truebndhighvec = ROOT.TVectorD(len(q2truebndhighlist))
        for i in range(len(q2truebndhighlist)):
            q2truebndhighvec[i]=float(q2truebndhighlist[i])
        xsec_norm_q2_true_bnd_high.AddAtAndExpand(q2truebndhighvec, nelem)


    nelem+=1 #Only have this line once par par

# Resize the vectors to a reasonable number of elements
xsec_param_nom_unnorm.ResizeTo(nelem)
xsec_param_nom.ResizeTo(nelem)
xsec_param_prior_unnorm.ResizeTo(nelem)
xsec_param_prior.ResizeTo(nelem)
xsec_param_lb.ResizeTo(nelem)
xsec_param_ub.ResizeTo(nelem)
xsec_stepscale.ResizeTo(nelem)
xsec_param_id.ResizeTo(nelem, 2)

xsec_cov = ROOT.TMatrixD(nelem, nelem)

# Loop over the elements
for i in range(nelem):

    # The diagonal entries in the covariance is simply the error squared
    xsec_cov[i][i] = xsec_error[i]*xsec_error[i]

    for item in corr_dicts[i]:
        # Script will fail if not in list of parameters; this is desired behavior. Could do more neatly, but meh
        index = corr_names.index(item)
        corr1 = corr_dicts[i][item]
        corr2 = 0;
        if (str(corr_names[i]) in corr_dicts[index].keys()):
            corr2 = corr_dicts[index][str(corr_names[i])]
        if (abs(corr1 - corr2) <= 1e-07):
            xsec_cov[i][index] = xsec_cov[index][i] = corr1 * xsec_error[i] * xsec_error[index]
        elif (corr2==0):
            print "Warning: no reciprocal correlation between "+str(item)+" and "+str(corr_names[i])+" proceeding with non-reciprocal"
            xsec_cov[i][index]=xsec_cov[index][i] = corr1 * xsec_error[i] * xsec_error[index]
        else:
            print "Correlations between "+str(item)+" and "+str(corr_names[i])+" don't match. Quitting!"
            print "corr 1 = ", corr1, "corr 2 = ", corr2
            quit()
        

# Make the TH2D covariance matrix
hcov = ROOT.TH2D("hcov", "", nelem, 0, nelem, nelem, 0, nelem)
# Set the content
for i in range(nelem):
    for j in range(nelem):
        if (xsec_cov[i][j] >= 0):
            hcov.SetBinContent(i+1, j+1, math.sqrt(xsec_cov[i][j]))
        else:
            hcov.SetBinContent(i+1, j+1, -1*math.sqrt(abs(xsec_cov[i][j])))

for i in range(nelem):
    hcov.GetXaxis().SetBinLabel(i+1, str(corr_names[i]))
    hcov.GetYaxis().SetBinLabel(i+1, str(corr_names[i]))

# Set the covariance matrix pretty
hcov.SetMaximum(1)
hcov.SetMinimum(-1)
hcov.GetXaxis().LabelsOption("v")
hcov.GetXaxis().SetLabelSize(0.03)
hcov.GetYaxis().SetLabelSize(0.03)
hcov.GetZaxis().SetTitle("#sqrt{|V_{ij}|}#times sign(V_{ij})")

# Write the output file
file = ROOT.TFile(sys.argv[2], "RECREATE")
xsec_param_names.Write("xsec_param_names", 1)
xsec_sk_spline_names.Write("sk_spline_names", 1);
xsec_sk_spline_modes.Write("sk_spline_modes",1);
xsec_nd_spline_names.Write("nd_spline_names", 1);
xsec_param_nom.Write("xsec_param_nom")
xsec_param_nom_unnorm.Write("xsec_param_nom_unnorm")
xsec_param_prior.Write("xsec_param_prior")
xsec_param_prior_unnorm.Write("xsec_param_prior_unnorm")
xsec_param_lb.Write("xsec_param_lb")
xsec_param_ub.Write("xsec_param_ub")
xsec_param_id.Write("xsec_param_id")
xsec_stepscale.Write("xsec_stepscale")
xsec_norm_modes.Write("xsec_norm_modes", 1)
xsec_norm_horncurrents.Write("xsec_norm_horncurrents", 1)
xsec_norm_etru_bnd_low.Write("xsec_norm_etru_bnd_low",1)
xsec_norm_etru_bnd_high.Write("xsec_norm_etru_bnd_high",1)
xsec_norm_q2_true_bnd_low.Write("xsec_norm_q2_true_bnd_low",1)
xsec_norm_q2_true_bnd_high.Write("xsec_norm_q2_true_bnd_high",1)
xsec_norm_elements.Write("xsec_norm_elements", 1)
xsec_norm_nupdg.Write("xsec_norm_nupdg", 1)
xsec_norm_prod_nupdg.Write("xsec_norm_prod_nupdg", 1)
xsec_cov.Write("xsec_cov")
hcov.Write("hcov")
file.Close()

#filename=sys.argv[1];
#filename = filename.replace('.xml', '.tex')
#print filename
#texfile = open(filename, 'w')

#for i in range(nelem):
#    snomu =  "%.2f" % xsec_param_nom_unnorm[i]
#    snom =  "%.2f" % xsec_param_nom[i]
#    sprioru =  "%.2f" % xsec_param_prior_unnorm[i]
#    sprior =  "%.2f" % xsec_param_prior[i]
#    serroru = ""
#    if(snomu==snom):
#        serroru =  "%.2f" % xsec_error[i];
#    else:    
#        serroru =  "%.2f" % (xsec_error[i]*xsec_param_nom[i])
#    serror =  "%.2f" % xsec_error[i]
#    slb =  "%.2f" % xsec_param_lb[i]
#    sub =  "%.2f" % xsec_param_ub[i]
#    texfile.write(sprioru+"&"+serroru+"&"+snomu+"&"+sprior+"&"+serror+"&"+slb+"&"+sub+"\\\\\n")

