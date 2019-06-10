#!/usr/bin/env python
import sys, getopt
import array
import re
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex

import skipper_utils

gROOT.SetBatch(True)
gStyle.SetOptStat(110011)
gStyle.SetOptFit(1)
gStyle.SetPalette(57)

EXPOSURE = 60.0 #minutes, per file

READOUT = 1.0*60+27 #minutes, per file

NROWS = 700 #number of rows read out, per file

OHDU = 1

gain = 300

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

infiles = remainder[1:]
infiles.sort(key=skipper_utils.decodeRunnum)

print(EXPOSURE,READOUT)
#if (len(sys.argv)<3):
#    print("not enough input args")
#    sys.exit(1)
outfilename = remainder[0].rpartition(".")[0] #strip off file extension, if any
if outfilename=="":
    outfilename = remainder[0]

data = TChain("skPixTree")
osm = TChain("osMeanTree")
numRuns = 0
for i in range(0,len(infiles)):
    print(infiles[i])
    data.Add(infiles[i])
    osm.Add(infiles[i])
    numRuns += 1
    thefile = TFile(infiles[i])
    header = thefile.Get("headerTree_0")
    ogl = skipper_utils.getHeaderValue(header,"OGAL")
    swl = skipper_utils.getHeaderValue(header,"SWAL")
    print("OGL={0}, SWL={1}".format(ogl,swl))


#x=0 and y=0 look weird (negative or extreme values)
#x=[370,450] is overscan
#x=[1,7] is prescan
#x=[8,369] looks real - matches dimension in paper (362 pixels)
#y=[625,700] is vertical overscan
#y=[1,624] looks real - matches dimension in paper (624 pixels)

funcformulas = []
for ipeak in range(0,5):
    funcformulas.append("[0]*(TMath::Gaus(x,[1]+{0}*[2],[2]*[3],1)*TMath::Poisson({0},[4]) )".format(ipeak))

fitfunc = skipper_utils.poissonFitfunc()

isMonsoon = (gain<100)

pixmax = 4.0*gain
pixmin = -1.0*gain
pixnbin = 50
pixbinning = "({0},{1},{2})".format(pixnbin,pixmin,pixmax)
binwidth = (pixmax-pixmin)/pixnbin

c = TCanvas("c","c",1200,900);
c.Print(outfilename+".pdf[")

c.SetLogy(1)



ohdus = []
ohduelists = []
fitvals = []
badFitList = []

for ohdu in range(1,5):
    elistname = "e{0}".format(ohdu)
    n = data.Draw(">>"+elistname,"ohdu=={0}".format(ohdu))
    if n>0:
        ohdus.append(ohdu)
        ohduelists.append(gDirectory.Get(elistname))

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

hists = []
c.Clear()
c.Divide(len(ohdus),1)
for iHdu in range(0,len(ohdus)):
    c.cd(iHdu+1)
    gPad.SetLogy(1)

    hname = "hdata{0}".format(ohdus[iHdu])

    data.SetEventList(ohduelists[iHdu])
    data.Draw("pix>>"+hname+pixbinning,"x>8 && y>0","colz")
    hists.append(gDirectory.Get(hname))

    fitvals.append({})
    badFitList.append(skipper_utils.fitPeaksGaus(hists[iHdu],gain,fitvals[iHdu]))

c.cd()
c.Print(outfilename+".pdf");

print(fitvals)

for iHdu in range(0,len(ohdus)):
    c.cd(iHdu+1)
    gPad.SetLogy(1)
    if not badFitList[iHdu]:
        s = skipper_utils.fitPeaksPoisson(fitfunc,hists[iHdu],fitvals[iHdu])
        badFitList[iHdu] = (int(s)!=0)

c.cd()
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

htypes = ["prescan","overscan_x","overscan_y","active"]
regioncuts = ["x>0 && x<8 && y>0","x>369 && y>0","x>=8 && x<=369 && y>624","x>=8 && x<=369 && y>0 && y<=624"]
exposures = [0.0, (362.0/450)*READOUT/NROWS, READOUT*(624.0/NROWS), 0.5*READOUT*(624.0/NROWS)+EXPOSURE]

latex = TLatex()
latex.SetNDC(True)
c.Clear()
c.Divide(len(ohdus),4)
histList = [] #not used, but needed so Python doesn't delete the histograms
for iHdu in range(0,len(ohdus)):
    data.SetEventList(ohduelists[iHdu])
    for hnum in range(0,4):
        c.cd(iHdu + hnum*len(ohdus) +1)
        gPad.SetLogy(1)
        hname = "h{0}{1}".format(ohdus[iHdu],htypes[hnum])
        data.Draw("pix>>"+hname+pixbinning,regioncuts[hnum],"colz")
        h = gDirectory.Get(hname)
        histList.append(h)
        if not badFitList[iHdu]:
            s = skipper_utils.fitPeaksPoisson(fitfunc,h,fitvals[iHdu])
            if (int(s)==0): #good fit
                #arrPoisson.append(s.Parameter(4))
                #arrPoissonErr.append(s.Error(4))
                #arrActiveTime.append(exposures[hnum])
                #arrZero.append(0.0)
                latex.DrawLatex(0.3,0.8,"mu={0} \pm {1}".format(s.Parameter(4),s.Error(4)))


c.cd()
c.Print(outfilename+".pdf");
c.Clear()

#if not badFit:
    #print(arrPoisson)
    #print(arrPoissonErr)
    #print(arrActiveTime)
    #c.SetLogy(0)
    #gStyle.SetStatX(0.4)
    #thegraph = TGraphErrors(len(arrPoisson),arrActiveTime,arrPoisson,arrZero,arrPoissonErr)
    #thegraph.Draw("AP")
    #thegraph.Fit("pol1")
    #thegraph.SetMinimum(0.0)
    #thegraph.SetMarkerSize(10)
    #c.Print(outfilename+".pdf");
#c.SetLogy(1)
#gStyle.SetStatX(0.9)

pixcuts = [("any","ele>0.6"),
        ("small","ele>0.6 && ele<3.5"),
        ("big","ele>10"),
        ("neg","ele<-0.7")]

gStyle.SetOptStat(0)
c.Clear()
c.Divide(len(ohdus),4)
histList = [] #not used, but needed so Python doesn't delete the histograms
for iHdu in range(0,len(ohdus)):
    data.SetEventList(ohduelists[iHdu])
    if not badFitList[iHdu]:
        data.SetAlias("ele","(pix-{0})/{1}".format(fitvals[iHdu]["zero"],fitvals[iHdu]["gain"]))

        #gPad.SetLogy(1)
        #data.Draw("ele>>hdata(500,-1,4)","x>8 && y>0 && "+ohducut,"colz")
        #c.Print(outfilename+".pdf");

        for hnum in range(0,4):
            c.cd(iHdu + hnum*len(ohdus) +1)
            hname = "h2d{0}_{1}".format(ohdus[iHdu],pixcuts[hnum][0])
            binning = "(450,-0.5,449.5,700,-0.5,699.5)"
            print data.Draw("y:x>>"+hname+binning,"x>0 && y>0 && "+pixcuts[hnum][1],"")
c.cd()
c.Print(outfilename+".pdf");


c.SetLogy(1)
gStyle.SetOptStat(0)

linfunc = TF1("linfunc","[0]*(x+[1])",1,624)
linfunc.FixParameter(1,(638.0*READOUT/EXPOSURE - 1))
#p0+p1*1 = x*READOUT
#p0+p1*639 = x*(EXPOSURE+READOUT)
#p1*638 = x*EXPOSURE
#p1 = x*EXPOSURE/638
#p0 = x*(READOUT - EXPOSURE/638)
#p0/p1 = 638*READOUT/EXPOSURE - 1
hists = []
c.Clear()
c.Divide(len(ohdus),3)
for iHdu in range(0,len(ohdus)):
    data.SetEventList(ohduelists[iHdu])
    #data.Draw("ele>>hdata(75,-2.5,22.5)","x>0 && y>0 && "+ohducut,"colz")
    #c.Print(outfilename+".pdf");

    data.Draw(">>elist","x>0 && y>0 && ele>0.6 && ele<3.5")
    data.SetEventList(gDirectory.Get("elist"))

    #data.Draw("ele>>hdata(500,0,4)","","colz")
    #c.Print(outfilename+".pdf");
    #c.SetLogy(0)

    c.cd(iHdu + 0*len(ohdus) +1)
    data.Draw("x>>h{0}x(450,-0.5,449.5)".format(ohdus[iHdu]),"","colz")

    #data.Draw("x>>hdata(45,-0.5,449.5)","","colz")
    #c.Print(outfilename+".pdf");

    #data.Draw("x>>hdata(200,359.5,379.5)","x>0 && y>0 && pix>250","colz")
    #c.Print(outfilename+".pdf");

    c.cd(iHdu + 1*len(ohdus) +1)
    data.Draw("y>>h{0}y(700,-0.5,699.5)".format(ohdus[iHdu]),"x>=8 && x<=369","colz")

    c.cd(iHdu + 2*len(ohdus) +1)
    data.Draw("y>>h{0}y2(26,0.5,624.5)".format(ohdus[iHdu]),"x>=8 && x<=369","colz")
    h = gDirectory.Get("h{0}y2".format(ohdus[iHdu]))

    h.Sumw2()
    h.Fit("linfunc","Q")
    gStyle.SetStatX(0.4)

c.cd()
c.Print(outfilename+".pdf");


    #data.Draw("y>>hdata(200,600,650)","x>=8 && x<=369","colz")
    #c.Print(outfilename+".pdf");

    #data.Draw("y>>hdata(70,-0.5,699.5)","x>=8 && x<=100","colz")
    #c.Print(outfilename+".pdf");

    #data.Draw("y>>hdata(70,-0.5,699.5)","x>369","colz")
    #c.Print(outfilename+".pdf");

c.Print(outfilename+".pdf]");
#c.Print(outfilename+".pdf]");
#sys.exit(0)
