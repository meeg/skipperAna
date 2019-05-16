#!/usr/bin/env python
import sys
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend

gROOT.SetBatch(True)
gStyle.SetOptStat(110011)
gStyle.SetOptFit(1)

if (len(sys.argv)<3):
    print("not enough input args")
    sys.exit(1)
outfilename = sys.argv[1]

data = TChain("skPixTree")
numRuns = 0
for i in range(2,len(sys.argv)):
    print(sys.argv[i])
    data.Add(sys.argv[i])
    numRuns += 1


#x=0 and y=0 look weird (negative or extreme values)
#x=[370,450] is overscan
#x=[1,7] has values similar to overscan
#x=[8,369] looks real - matches dimension in paper (362 pixels)
#y=[625,700] is overscan
#y=[1,624] looks real - matches dimension in paper (624 pixels)

#data.Add("skp_1hr_33_12.root")
#data.Add("skp_1hr_34_12.root")
#data.Add("skp_1hr_35_12.root")
#data.Add("skp_1hr_36_12.root")
#data.Add("skp_1hr_37_12.root")

#outfilename="backgroundfit_histos"

EXPOSURE = 1.0/24 #days, per file

READOUT = (2*60+13)/60/24 #days, per file

#dark current for real pixel: exposure+readout
#for overscan: 
ipeak = 0
#peakfuncs = []
#funcnames = []
funcformulas = []
for ipeak in range(0,5):
    funcformulas.append("[0]*(TMath::Gaus(x,[1]+{0}*[2],[2]*[3],1)*TMath::Poisson({0},[4]) )".format(ipeak))
    #peakfuncs.append(TF1("peak{0}".format(ipeak),"[0]*(TMath::Gaus(x,[1]+{0}*[2],[2]*[3])*TMath::Poisson({0},[4]) )".format(ipeak)))
    #funcnames.append("peak{0}".format(ipeak))

#print(" + ".join(funcnames))
fitfunc = TF1("fitfunc"," + ".join(funcformulas))

#fitfunc = TF1("fitfunc","[0]*(TMath::Gaus(x,[1],[3])*TMath::Poisson(0,[4]) + TMath::Gaus(x,[1]+[2],[3])*TMath::Poisson(1,[4]) + TMath::Gaus(x,[1]+2*[2],[3])*TMath::Poisson(2,[4]) )",-2,3)


c = TCanvas("c","c",1200,900);
c.Print(outfilename+".pdf[")

c.SetLogy(1)
data.Draw("pix>>hdata(500,-500,2000)","x>0 && y>0","colz")
hdata = gDirectory.Get("hdata")
fitfunc.SetParameters(hdata.GetEntries()*5,0,500,0.1,0.1)
s = hdata.Fit(fitfunc,"S","")
s = hdata.Fit(fitfunc,"SL","")
c.Print(outfilename+".pdf");

c.Clear()
c.Divide(1,2)
c.cd(1)
gPad.SetLogy(1)
data.Draw("pix>>hoverscan(500,-500,2000)","x>369 && y>0","colz")
hoverscan = gDirectory.Get("hoverscan")
fitfunc.SetParameters(hoverscan.GetEntries()*5,0,500,0.1,0.01)
s = hoverscan.Fit(fitfunc,"S","")
s = hoverscan.Fit(fitfunc,"SL","")
#c.Print(outfilename+".pdf");

c.cd(2)
gPad.SetLogy(1)
data.Draw("pix>>hactive(500,-500,2000)","x>=8 && x<=369 && y>0 && y<=624","colz")
hactive = gDirectory.Get("hactive")
fitfunc.SetParameters(hactive.GetEntries()*5,0,500,0.1,0.1)
s = hactive.Fit(fitfunc,"S","")
s = hactive.Fit(fitfunc,"SL","")
c.cd()
c.Print(outfilename+".pdf");
c.Clear()

s0 = hdata.Fit("gaus","S","",-200,200)
#c.Print(outfilename+".pdf");

s1 = hdata.Fit("gaus","S","",300,700)
#c.Print(outfilename+".pdf");

zero = s0.Parameter(1)
aduPerElectron = s1.Parameter(1) - s0.Parameter(1)

s0 = hdata.Fit("gaus","S","",zero-0.4*aduPerElectron,zero+0.4*aduPerElectron)
c.Print(outfilename+".pdf");

s1 = hdata.Fit("gaus","S","",zero+0.7*aduPerElectron,zero+1.4*aduPerElectron)
c.Print(outfilename+".pdf");

zero = s0.Parameter(1)
aduPerElectron = s1.Parameter(1) - s0.Parameter(1)

data.SetAlias("ele","(pix-{0})/{1}".format(zero,aduPerElectron))

data.Draw("ele>>hdata(500,-1,4)","x>0 && y>0","colz")
c.Print(outfilename+".pdf");

c.SetLogy(0)
data.Draw("y:x>>h2d(450,-0.5,449.5,700,-0.5,699.5)","x>0 && y>0 && ele>3.5","colz")
c.Print(outfilename+".pdf");

c.SetLogy(1)
data.Draw("ele>>hdata(75,-2.5,22.5)","x>0 && y>0","colz")
c.Print(outfilename+".pdf");

data.Draw(">>elist","x>0 && y>0 && ele>0.6 && ele<3.5")
data.SetEventList(gDirectory.Get("elist"))

data.Draw("ele>>hdata(500,0,4)","","colz")
c.Print(outfilename+".pdf");
c.SetLogy(0)

data.Draw("x>>hdata(450,-0.5,449.5)","","colz")
c.Print(outfilename+".pdf");

data.Draw("x>>hdata(45,-0.5,449.5)","","colz")
c.Print(outfilename+".pdf");

#data.Draw("x>>hdata(200,359.5,379.5)","x>0 && y>0 && pix>250","colz")
#c.Print(outfilename+".pdf");

data.Draw("y>>hdata(70,-0.5,699.5)","x>=8 && x<=369","colz")
c.Print(outfilename+".pdf");

#data.Draw("y>>hdata(200,600,650)","x>=8 && x<=369","colz")
#c.Print(outfilename+".pdf");

data.Draw("y>>hdata(70,-0.5,699.5)","x>=8 && x<=100","colz")
c.Print(outfilename+".pdf");

data.Draw("y>>hdata(70,-0.5,699.5)","x>369","colz")
c.Print(outfilename+".pdf");

c.Print(outfilename+".pdf]");
