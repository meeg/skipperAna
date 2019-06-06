import sys, getopt
import array
import re
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex

def getHeaderValue(tree,name):
    tree.GetEntry(0)
    stringdata = tree.GetBranch(name).GetLeaf("string").GetValuePointer()
    b = bytearray()
    for iword in range(0,len(stringdata)):
        word = stringdata[iword]
        #print(word)
        for ibyte in range(0,8):
            b.append(word & 0xFF)
            word >>= 8
    #print(b)
    return float(b)

def decodeRunnum(filename):
    ints = [int(s) for s in filename.split("/")[-1].rpartition(".")[0].split("_") if s.isdigit()] #get the basename, strip off the file extension, extract underscore-separated ints
    return ints[0]

def formatVoltage(val):
    if (val<0):
        val *= -1
    if (round(val)==val):
        val = int(val)
    return str(val).replace(".","p")

def fitPeaksGaus(h,gain,fitvals):
    s0 = h.Fit("gaus","QS","",-0.4*gain,0.4*gain)
    s1 = h.Fit("gaus","QS","",0.7*gain,1.4*gain)
    badFit = (int(s0)!=0 or int(s1)!=0)
    if not badFit:
        fitvals["zero"] = s0.Parameter(1)
        fitvals["gain"] = s1.Parameter(1) - s0.Parameter(1)
        fitvals["noise"] = s0.Parameter(2)/fitvals["gain"]
        badFit = (fitvals["gain"]<0.5*gain or fitvals["gain"]>2*gain or fitvals["noise"]<0.01 or fitvals["noise"]>1)
    if not badFit:
        s0 = h.Fit("gaus","QS","",fitvals["zero"]-0.4*fitvals["gain"],fitvals["zero"]+0.4*fitvals["gain"])
        s1 = h.Fit("gaus","QS+","",fitvals["zero"]+0.7*fitvals["gain"],fitvals["zero"]+1.4*fitvals["gain"])
        badFit = (int(s0)!=0 or int(s1)!=0)
    if not badFit:
        fitvals["zero"] = s0.Parameter(1)
        fitvals["gain"] = s1.Parameter(1) - s0.Parameter(1)
        fitvals["noise"] = s0.Parameter(2)/fitvals["gain"]
        print("zero={0}, gain={1}, noise={2}".format(fitvals["zero"],fitvals["gain"],fitvals["noise"]))
        badFit = (fitvals["gain"]<0.5*gain or fitvals["gain"]>2*gain or fitvals["noise"]<0.01 or fitvals["noise"]>1)
    return badFit

def fitPeaksPoisson(fitfunc,h,fitvals):
    fitfunc.SetParameters(h.GetEntries()*h.GetXaxis().GetBinWidth(1),fitvals["zero"],fitvals["gain"],fitvals["noise"],0.01)
    h.Fit(fitfunc,"QS","",fitvals["zero"]-0.5*fitvals["gain"],fitvals["zero"]+2.5*fitvals["gain"])
    return h.Fit(fitfunc,"QSL","",fitvals["zero"]-0.5*fitvals["gain"],fitvals["zero"]+2.5*fitvals["gain"])

def poissonFitfunc():
    funcformulas = []
    for ipeak in range(0,5):
        funcformulas.append("[0]*(TMath::Gaus(x,[1]+{0}*[2],[2]*[3],1)*TMath::Poisson({0},[4]) )".format(ipeak))

    return TF1("poisson_fitfunc"," + ".join(funcformulas))
