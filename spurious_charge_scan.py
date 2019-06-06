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

OHDU = 1

gain = 400

settingNames=["SWAL","OGAL"]

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'g:o:x:y:h')

for opt, arg in options:
    if opt=='-g':
        gain = float(arg)
    if opt=='-o':
        OHDU = int(arg)
    if opt=='-x':
        settingNames[0] = arg
    if opt=='-y':
        settingNames[1] = arg
    elif opt=='-h':
        print("\nUsage: "+sys.argv[0]+" <output basename> <root files>")
        print("Arguments: ")
        print("\t-g, --gain: initial guess for gain (ADU per electron)")
        print("\t-o, --ohdu: HDU to analyze")
        print("\n")
        sys.exit(0)


if (len(remainder)<2):
    print(sys.argv[0]+' <output basename> <root files>')
    sys.exit()

outfilename = remainder[0]
infiles = remainder[1:]
infiles.sort(key=skipper_utils.decodeRunnum)

fitfunc = skipper_utils.poissonFitfunc()

isMonsoon = (gain<100)

pixmax = 4.0*gain
pixmin = -1.0*gain
pixnbin = 50
pixbinning = "({0},{1},{2})".format(pixnbin,pixmin,pixmax)
binwidth = (pixmax-pixmin)/pixnbin

c = TCanvas("c","c",1200,900);
c.Print(outfilename+".pdf[")

outfile = TFile(outfilename+".root","RECREATE")

latex = TLatex()
latex.SetNDC(True)

#x=0 and y=0 look weird (negative or extreme values)
#x=[370,450] is overscan
#x=[1,7] is prescan
#x=[8,369] looks real - matches dimension in paper (362 pixels)
#y=[625,700] is vertical overscan
#y=[1,624] looks real - matches dimension in paper (624 pixels)

htypes = ["allpix","prescan","overscan","active"]
regioncuts = ["x>0 && y>0","x>0 && x<8 && y>0","x>369 && y>0","x>=8 && x<=369 && y>0 && y<=624"]

settingListDict={}
for settingName in settingNames:
    settingListDict[settingName]=[]
numRuns = 0
settingString=""#compact string for histogram names
settingPrintString=""#readable string for labels
histDict={}#setting string to list of list of histograms
dataChain = TChain("skPixTree")
for iFile in range(0,len(infiles)):
    prevSettingString = settingString
    print(infiles[iFile])
    runnum = skipper_utils.decodeRunnum(infiles[iFile])
    thefile = TFile(infiles[iFile])
    header = thefile.Get("headerTree_0")
    data = thefile.Get("skPixTree")
    outfile.cd()

    settingString=""
    settingPrintList=[]
    for settingName in settingNames:
        settingVal = skipper_utils.getHeaderValue(header,settingName)
        settingPrintList.append("{0}={1}".format(settingName,settingVal))
        settingWord = "{0}{1}".format(settingName,skipper_utils.formatVoltage(settingVal))
        settingString += settingWord
        if settingListDict[settingName].count(settingWord)==0:
            settingListDict[settingName].append(settingWord)
    settingPrintString = ", ".join(settingPrintList)
    print(settingPrintString)

    if settingString not in histDict:
        histDict[settingString]=[]
    histsThisFile = []
    histDict[settingString].append(histsThisFile)
    fileHistString = settingString+"FILE{0}".format(runnum)

    c.Clear()
    c.Divide(1,4)
    latex.DrawLatex(0.1,0.9,settingPrintString)

    fitvals={}

    for hnum in range(0,4):
        c.cd(hnum+1)
        gPad.SetLogy(1)
        hname = "h"+htypes[hnum]+fileHistString
        data.Draw("pix>>"+hname+pixbinning,regioncuts[hnum],"colz")
        h = gDirectory.Get(hname)
        histsThisFile.append(h)
        if hnum==0:
            badFit = skipper_utils.fitPeaksGaus(h,gain,fitvals)
        if not badFit:
            skipper_utils.fitPeaksPoisson(fitfunc,h,fitvals)

    c.cd()
    c.Print(outfilename+".pdf");

sumHistDict={}#setting string to list of histograms
gStyle.SetOptStat(0)
c.Clear()
c.Divide(len(settingListDict[settingNames[0]]),len(settingListDict[settingNames[1]]))
for iHist in range(0,4):
    settingStringList=["" for x in settingNames]
    for settingWord0 in settingListDict[settingNames[0]]:
        for settingWord1 in settingListDict[settingNames[1]]:
            settingString = settingWord0+settingWord1
            #print(settingString)
            if settingString not in sumHistDict:
                sumHistDict[settingString]=[]
            #print(histDict[settingString])
            sumHist = histDict[settingString][0][iHist].Clone(htypes[iHist]+settingString)
            sumHistDict[settingString].append(sumHist)
            sumHist.Reset()
            for fileHists in histDict[settingString]:
                sumHist.Add(fileHists[iHist])
            canvasNum = settingListDict[settingNames[0]].index(settingWord0) + settingListDict[settingNames[1]].index(settingWord1)*len(settingListDict[settingNames[0]]) + 1
            #print(canvasNum)
            c.cd(canvasNum)
            gPad.SetLogy(1)
            sumHist.Draw()
            fitvals={}
            badFit = skipper_utils.fitPeaksGaus(sumHist,gain,fitvals)
            if not badFit:
                skipper_utils.fitPeaksPoisson(fitfunc,sumHist,fitvals)
            latex.DrawLatex(0.3,0.8,settingString)
    c.cd()
    c.Print(outfilename+".pdf");

print settingListDict
c.Print(outfilename+".pdf]");
outfile.Write()
outfile.Close()
