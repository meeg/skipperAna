#!/usr/bin/env python
import sys, getopt
import array
import re
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex

gROOT.SetBatch(True)
gStyle.SetOptStat(110011)
gStyle.SetOptFit(1)
gStyle.SetPalette(57)

EXPOSURE = 60.0 #minutes, per file

READOUT = 1.0*60+27 #minutes, per file

NROWS = 700 #number of rows read out, per file

OHDU = 1

gain = 500

OHDPLOTS = False
OSMPLOTS = False

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'e:r:g:o:h', ['exposure','readout','gain','ohdu','help'])

for opt, arg in options:
    if opt in ('-e', '--exposure'):
        EXPOSURE = float(arg)
    if opt in ('-r', '--readout'):
        READOUT = float(arg)
    if opt in ('-g', '--gain'):
        gain = float(arg)
    if opt in ('-o', '--ohdu'):
        OHDU = int(arg)
    elif opt in ('-h', '--help'):
        print "\nUsage: "+sys.argv[0]+" <output basename> <root files>"
        print "Arguments: "
        print "\t-e, --exposure: exposure time (minutes) for normalization"
        print "\t-r, --readout: readout time (minutes) for normalization"
        print "\t-g, --gain: initial guess for gain (ADU per electron)"
        print "\t-o, --ohdu: HDU to analyze"
        print "\n"
        sys.exit(0)


if (len(remainder)<2):
    print sys.argv[0]+' <output basename> <root files>'
    sys.exit()

def runnum(filename):
    ints = [int(c) for c in re.split(r'(\d+)',filename) if c.isdigit()]
    return ints[-2]

remainder.sort(key=runnum)
#for i in range(1,len(remainder)):
    #print([int(c) for c in re.split(r'(\d+)',remainder[i]) if c.isdigit()])
    #print runnum(remainder[i])


print(EXPOSURE,READOUT)
#if (len(sys.argv)<3):
#    print("not enough input args")
#    sys.exit(1)
outfilename = remainder[0]

data = TChain("skPixTree")
osm = TChain("osMeanTree")
numRuns = 0
for i in range(1,len(remainder)):
    print(remainder[i])
    data.Add(remainder[i])
    osm.Add(remainder[i])
    numRuns += 1

#x=0 and y=0 look weird (negative or extreme values)
#x=[370,450] is overscan
#x=[1,7] is prescan
#x=[8,369] looks real - matches dimension in paper (362 pixels)
#y=[625,700] is vertical overscan
#y=[1,624] looks real - matches dimension in paper (624 pixels)

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

isMonsoon = (gain<100)

pixmax = 4.0*gain
pixmin = -1.0*gain
pixnbin = 250
pixbinning = "({0},{1},{2})".format(pixnbin,pixmin,pixmax)
binwidth = (pixmax-pixmin)/pixnbin

c = TCanvas("c","c",1200,900);
c.Print(outfilename+".pdf[")

c.SetLogy(1)

if (OSMPLOTS):
    if isMonsoon:
        osm.Draw("osm","spl>20")
    else:
        osm.Draw("osm","spl>5")
    c.Print(outfilename+".pdf");

if (OHDPLOTS):
    c.Clear()
    c.Divide(1,5)
    for i in range(1,6):
        c.cd(i)
        gPad.SetLogy(1)
        data.Draw("pix>>h{0}{1}".format(i,pixbinning),"x>8 && y>0 && ohdu=={0}".format(i),"colz")
        h = gDirectory.Get("h{0}".format(i))
    c.cd()
    c.Print(outfilename+".pdf");
    for i in range(1,6):
        c.cd(i)
        gPad.SetLogy(1)
        data.Draw("pix>>h{0}{1}".format(i,pixbinning),"x>369 && y>0 && ohdu=={0}".format(i),"colz")
        h = gDirectory.Get("h{0}".format(i))
    c.cd()
    c.Print(outfilename+".pdf");
    c.Clear()

ohducut = "ohdu=={0}".format(OHDU)

data.Draw("pix>>hdata"+pixbinning,"x>8 && y>0 && "+ohducut,"colz")
hdata = gDirectory.Get("hdata")

#c.Print(outfilename+".pdf");
#c.Print(outfilename+".pdf]");
#sys.exit(0)

s0 = hdata.Fit("gaus","QS","",-0.4*gain,0.4*gain)
s1 = hdata.Fit("gaus","QS","",0.6*gain,1.4*gain)
aduZero = s0.Parameter(1)
aduPerElectron = s1.Parameter(1) - s0.Parameter(1)

s0 = hdata.Fit("gaus","QS","",aduZero-0.4*aduPerElectron,aduZero+0.4*aduPerElectron)
c.Print(outfilename+".pdf");
s1 = hdata.Fit("gaus","QS","",aduZero+0.7*aduPerElectron,aduZero+1.4*aduPerElectron)
c.Print(outfilename+".pdf");
aduZero = s0.Parameter(1)
aduPerElectron = s1.Parameter(1) - s0.Parameter(1)
noise = s0.Parameter(2)/aduPerElectron




fitfunc.SetParameters(hdata.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.1)
s = hdata.Fit(fitfunc,"QS","")
s = hdata.Fit(fitfunc,"QSL","")
c.Print(outfilename+".pdf");

#data.Draw("pix>>hdata"+pixbinning,"x>0 && y>0 && "+ohducut,"colz")
#hdata = gDirectory.Get("hdata")
#fitfunc.SetParameters(hdata.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.1)
#s = hdata.Fit(fitfunc,"QS","")
#s = hdata.Fit(fitfunc,"QSL","")
#c.Print(outfilename+".pdf");

arrPoisson = array.array('d')
arrPoissonErr = array.array('d')
arrActiveTime = array.array('d')
arrZero = array.array('d')
#arrExposure = array.array('d')
#arrReadout = array.array('d')
#arrShift = array.array('d')

latex = TLatex()
c.Clear()
c.Divide(1,4)

c.cd(1)
gPad.SetLogy(1)
data.Draw("pix>>hprescan"+pixbinning,"x>0 && x<8 && y>0 &&"+ohducut,"colz")
prescan = gDirectory.Get("hprescan")
fitfunc.SetParameters(prescan.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.01)
prescan.Fit(fitfunc,"QS","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
s = prescan.Fit(fitfunc,"QSL","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
arrPoisson.append(s.Parameter(4))
arrPoissonErr.append(s.Error(4))
arrActiveTime.append(0.0)
arrZero.append(0.0)
latex.DrawLatex(500,100,"mu={0} \pm {1}".format(s.Parameter(4),s.Error(4)))

c.cd(2)
gPad.SetLogy(1)
data.Draw("pix>>hoverscan_x"+pixbinning,"x>369 && y>0 &&"+ohducut,"colz")
hoverscan_x = gDirectory.Get("hoverscan_x")
fitfunc.SetParameters(hoverscan_x.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.01)
hoverscan_x.Fit(fitfunc,"QS","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
s = hoverscan_x.Fit(fitfunc,"QSL","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
arrPoisson.append(s.Parameter(4))
arrPoissonErr.append(s.Error(4))
arrActiveTime.append((362.0/450)*READOUT/NROWS)
arrZero.append(0.0)
latex.DrawLatex(500,100,"mu={0} \pm {1}".format(s.Parameter(4),s.Error(4)))

c.cd(3)
gPad.SetLogy(1)
data.Draw("pix>>hoverscan_y"+pixbinning,"x>=8 && x<=369 && y>624 &&"+ohducut,"colz")
hoverscan_y = gDirectory.Get("hoverscan_y")
fitfunc.SetParameters(hoverscan_y.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.01)
hoverscan_y.Fit(fitfunc,"QS","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
s = hoverscan_y.Fit(fitfunc,"QSL","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
arrPoisson.append(s.Parameter(4))
arrPoissonErr.append(s.Error(4))
arrActiveTime.append(READOUT*(624.0/NROWS))
arrZero.append(0.0)
latex.DrawLatex(500,100,"mu={0} \pm {1}".format(s.Parameter(4),s.Error(4)))

c.cd(4)
gPad.SetLogy(1)
data.Draw("pix>>hactive"+pixbinning,"x>=8 && x<=369 && y>0 && y<=624 &&"+ohducut,"colz")
hactive = gDirectory.Get("hactive")
fitfunc.SetParameters(hactive.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.1)
hactive.Fit(fitfunc,"QS","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
s = hactive.Fit(fitfunc,"QSL","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
arrPoisson.append(s.Parameter(4))
arrPoissonErr.append(s.Error(4))
arrActiveTime.append(0.5*READOUT*(624.0/NROWS)+EXPOSURE)
arrZero.append(0.0)
latex.DrawLatex(500,100,"mu={0} \pm {1}".format(s.Parameter(4),s.Error(4)))

#print(s_act.Parameter(4), s_osy.Parameter(4), s_osy.Parameter(4)/(s_act.Parameter(4)-s_osy.Parameter(4)))

c.cd()
c.Print(outfilename+".pdf");
c.Clear()

#print(0.5*READOUT*(624.0/NROWS)+EXPOSURE)
#print(0.5*READOUT*(624.0/NROWS)+EXPOSURE)

print(arrPoisson)
print(arrPoissonErr)
print(arrActiveTime)
c.SetLogy(0)
gStyle.SetStatX(0.4)
thegraph = TGraphErrors(len(arrPoisson),arrActiveTime,arrPoisson,arrZero,arrPoissonErr)
thegraph.Draw("AP")
thegraph.Fit("pol1")
thegraph.SetMinimum(0.0)
thegraph.SetMarkerSize(10)
c.Print(outfilename+".pdf");
c.SetLogy(1)
gStyle.SetStatX(0.9)


data.SetAlias("ele","(pix-{0})/{1}".format(aduZero,aduPerElectron))

data.Draw("ele>>hdata(500,-1,4)","x>0 && y>0 && "+ohducut,"colz")
c.Print(outfilename+".pdf");

c.SetLogy(0)
gStyle.SetOptStat(0)
c.Clear()
c.Divide(2,2)
c.cd(1)
data.Draw("y:x>>h2d_any(450,-0.5,449.5,700,-0.5,699.5)","x>0 && y>0 && ele>0.6 && "+ohducut,"colz")

c.cd(2)
data.Draw("y:x>>h2d_small(450,-0.5,449.5,700,-0.5,699.5)","x>0 && y>0 && ele>0.6 && ele<3.5 && "+ohducut,"colz")

c.cd(3)
data.Draw("y:x>>h2d_big(450,-0.5,449.5,700,-0.5,699.5)","x>0 && y>0 && ele>3.5 && "+ohducut,"colz")

c.cd(4)
data.Draw("y:x>>h2d_neg(450,-0.5,449.5,700,-0.5,699.5)","x>0 && y>0 && ele<-0.7 && "+ohducut,"colz")

c.cd()
c.Print(outfilename+".pdf");
c.Clear()
c.SetLogy(1)
gStyle.SetOptStat(110011)
#data.Draw("ele>>hdata(75,-2.5,22.5)","x>0 && y>0 && "+ohducut,"colz")
#c.Print(outfilename+".pdf");

data.Draw(">>elist","x>0 && y>0 && ele>0.6 && ele<3.5 && "+ohducut)
data.SetEventList(gDirectory.Get("elist"))

data.Draw("ele>>hdata(500,0,4)","","colz")
c.Print(outfilename+".pdf");
c.SetLogy(0)

data.Draw("x>>hdata(450,-0.5,449.5)","","colz")
c.Print(outfilename+".pdf");

#data.Draw("x>>hdata(45,-0.5,449.5)","","colz")
#c.Print(outfilename+".pdf");

#data.Draw("x>>hdata(200,359.5,379.5)","x>0 && y>0 && pix>250","colz")
#c.Print(outfilename+".pdf");

data.Draw("y>>hdata(700,-0.5,699.5)","x>=8 && x<=369","colz")
c.Print(outfilename+".pdf");

data.Draw("y>>hdata(26,0.5,624.5)","x>=8 && x<=369","colz")
h = gDirectory.Get("hdata")

linfunc = TF1("linfunc","[0]*(x+[1])",1,624)
linfunc.FixParameter(1,(638.0*READOUT/EXPOSURE - 1))
#linfunc = TF1("linfunc","pol1",1,624)
#linfunc.SetParameters(h.Integral()/26/(1 + ((1+639)*0.5)/(638*READOUT/EXPOSURE - 1)),h.Integral()/26/(638*READOUT/EXPOSURE - 1 + (1+639)*0.5))
h.Sumw2()
h.Fit("linfunc")
#p0+p1*1 = x*READOUT
#p0+p1*639 = x*(EXPOSURE+READOUT)
#p1*638 = x*EXPOSURE
#p1 = x*EXPOSURE/638
#p0 = x*(READOUT - EXPOSURE/638)
#p0/p1 = 638*READOUT/EXPOSURE - 1
gStyle.SetStatX(0.4)

c.Print(outfilename+".pdf");

#data.Draw("y>>hdata(200,600,650)","x>=8 && x<=369","colz")
#c.Print(outfilename+".pdf");

#data.Draw("y>>hdata(70,-0.5,699.5)","x>=8 && x<=100","colz")
#c.Print(outfilename+".pdf");

#data.Draw("y>>hdata(70,-0.5,699.5)","x>369","colz")
#c.Print(outfilename+".pdf");

c.Print(outfilename+".pdf]");
