#!/usr/bin/python

import xml.etree.ElementTree as ET
import ROOT
import math
import sys

if len(sys.argv) != 3:
    print "Sorry, I need two arguments"
    print "./makeXSecMatrix.py input.xml output.root"
    sys.exit()

tree = ET.parse(sys.argv[1])
fromxml = tree.getroot()

maxelements=20

osc_param_names        = ROOT.TObjArray()
osc_param_nom          = ROOT.TVectorD(maxelements)
osc_error              = ROOT.TVectorD(maxelements)
osc_stepscale          = ROOT.TVectorD(maxelements)
osc_sigma              = ROOT.TVectorD(maxelements)
osc_flat_prior         = ROOT.TVectorD(maxelements)

osc_baseline = ROOT.TVectorD(1)
osc_density = ROOT.TVectorD(1)


nelem = 0
for child in fromxml:

    if( child.tag == 'parameter'):
	    # Get the name attribute
        name = ROOT.TObjString(child.attrib['name'])

	    # Add the name to the TObjArray
        osc_param_names.AddLast(name)
        osc_param_nom[nelem] = float(child.attrib['nom'])
        osc_error[nelem] = float(child.attrib['error'])
        osc_stepscale[nelem] = float(child.attrib['stepscale'])
        osc_sigma[nelem] = float(child.attrib['sigma'])
        osc_flat_prior[nelem] = float(child.attrib['FlatPrior'])

        nelem+=1 #Only have this line once par par

    if( child.tag == 'experiment'):
	    #KS: sligtlhy hardcoded but we can have only one baselin and earth density
        print "Setting baseline to be"
        print(float(child.attrib['L']))
        osc_baseline[0] = float(child.attrib['L'])
        osc_density[0] = float(child.attrib['density'])

# Resize the vectors to a reasonable number of elements
osc_param_nom.ResizeTo(nelem)
osc_error.ResizeTo(nelem)
osc_stepscale.ResizeTo(nelem)
osc_sigma.ResizeTo(nelem)
osc_flat_prior.ResizeTo(nelem)
osc_cov = ROOT.TMatrixD(nelem, nelem)

# Loop over the elements
for i in range(nelem):
	# The diagonal entries in the covariance is simply the error squared
	osc_cov[i][i] = osc_error[i]*osc_error[i]


# Make the TH2D covariance matrix
hcov = ROOT.TH2D("hcov", "", nelem, 0, nelem, nelem, 0, nelem)
# Set the content
for i in range(nelem):
    for j in range(nelem):
	    if (osc_cov[i][j] >= 0):
	        hcov.SetBinContent(i+1, j+1, math.sqrt(osc_cov[i][j]))
	    else:
	        hcov.SetBinContent(i+1, j+1, -1*math.sqrt(abs(osc_cov[i][j])))

for i in range(nelem):
    hcov.GetXaxis().SetBinLabel(i+1, str(osc_param_names[i]))
    hcov.GetYaxis().SetBinLabel(i+1, str(osc_param_names[i]))

# Set the covariance matrix pretty
hcov.SetMaximum(1)
hcov.SetMinimum(-1)
hcov.GetXaxis().LabelsOption("v")
hcov.GetXaxis().SetLabelSize(0.03)
hcov.GetYaxis().SetLabelSize(0.03)
hcov.GetZaxis().SetTitle("#sqrt{|V_{ij}|}#times sign(V_{ij})")

# Write the output file
file = ROOT.TFile(sys.argv[2], "RECREATE")
osc_param_names.Write("osc_param_names", 1)
osc_param_nom.Write("osc_nom")
osc_stepscale.Write("osc_stepscale")
osc_sigma.Write("osc_sigma")
osc_flat_prior.Write("osc_flat_prior")
osc_cov.Write("osc_cov")

osc_baseline.Write("osc_baseline")
osc_density.Write("osc_density")
hcov.Write("hcov")
file.Close()
