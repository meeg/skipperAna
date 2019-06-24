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

gain = 400

settingNames=["SWAL","OGAL"]

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'g:x:y:h')

for opt, arg in options:
    if opt=='-g':
        gain = float(arg)
    if opt=='-x':
        settingNames[0] = arg
    if opt=='-y':
        settingNames[1] = arg
    elif opt=='-h':
        print("\nUsage: "+sys.argv[0]+" <output basename> <root files>")
        print("Arguments: ")
        print("\t-g: initial guess for gain (ADU per electron)")
        print("\t-x: header variable to plot on the X axis")
        print("\t-y: header variable to plot on the Y axis")
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

ohdus = []

htypes = ["allpix","prescan","overscan","active"]
regioncuts = ["x>0 && y>0","x>0 && x<8 && y>0","x>369 && y>0","x>=8 && x<=369 && y>0 && y<=624"]

settingListDict={}
for settingName in settingNames:
    settingListDict[settingName]=[]

numRuns = 0
settingString=""#compact string for histogram names
histDict={}#setting string to list of list of histograms

for iFile in range(0,len(infiles)):
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
    histsThisFile = {}
    histDict[settingString].append(histsThisFile)
    fileHistString = settingString+"FILE{0}".format(runnum)

    ohduelists = []
    if (iFile==0): #first file: get the list of HDUs
        for ohdu in range(0,6):#LTA data has 1 through 4, Monsoon data has 2 through 5
            elistname = "e{0}".format(ohdu)
            n = data.Draw(">>"+elistname,"ohdu=={0}".format(ohdu))
            if n>0:
                ohdus.append(ohdu)
                ohduelists.append(gDirectory.Get(elistname))
    else: #use the list of HDUs, regenerate the event lists
        for iHdu in range(0,len(ohdus)):
            ohdu = ohdus[iHdu]
            elistname = "e{0}".format(ohdu)
            n = data.Draw(">>"+elistname,"ohdu=={0}".format(ohdu))
            ohduelists.append(gDirectory.Get(elistname))

    c.Clear()
    c.Divide(len(ohdus),len(htypes))
    latex.DrawLatex(0.1,0.95,settingPrintString)
    for iHdu in range(0,len(ohdus)):
        data.SetEventList(ohduelists[iHdu])
        fitvals={}
        for iHist in range(0,len(htypes)):
            c.cd(iHdu + iHist*len(ohdus) +1)
            #c.cd(iHist+1)
            gPad.SetLogy(1)
            hname = "h{0}{1}FILE{2}OHDU{3}".format(htypes[iHist],settingString,runnum,ohdus[iHdu])
            data.Draw("pix>>"+hname+pixbinning,regioncuts[iHist],"colz")
            h = gDirectory.Get(hname)
            histString = "h{0}OHDU{1}".format(htypes[iHist],ohdus[iHdu])
            histsThisFile[histString] = h
            if iHist==0:
                badFit = skipper_utils.fitPeaksGaus(h,gain,fitvals)
            if not badFit:
                skipper_utils.fitPeaksPoisson(fitfunc,h,fitvals)
        c.cd(iHdu + (len(htypes)-1)*len(ohdus) +1)
        latex.DrawLatex(0.4,0.025,"OHDU{0}".format(ohdus[iHdu]))
    c.cd()
    c.Print(outfilename+".pdf");

gStyle.SetOptStat(0)
for iHdu in range(0,len(ohdus)):
    for iHist in range(0,4):
        c.Clear()
        c.Divide(len(settingListDict[settingNames[0]]),len(settingListDict[settingNames[1]]))
        sumHistList=[]#list of histograms
        settingStringList=["" for x in settingNames]
        for settingWord0 in settingListDict[settingNames[0]]:
            for settingWord1 in settingListDict[settingNames[1]]:
                settingString = settingWord0+settingWord1
                hname = htypes[iHist]+settingString
                histString = "h{0}OHDU{1}".format(htypes[iHist],ohdus[iHdu])
                sumHist = histDict[settingString][0][histString].Clone(hname)
                sumHistList.append(sumHist)
                sumHist.Reset()
                for fileHists in histDict[settingString]:
                    sumHist.Add(fileHists[histString])
                canvasNum = settingListDict[settingNames[0]].index(settingWord0) + settingListDict[settingNames[1]].index(settingWord1)*len(settingListDict[settingNames[0]]) + 1
                c.cd(canvasNum)
                gPad.SetLogy(1)
                sumHist.Draw()
                fitvals={}
                badFit = skipper_utils.fitPeaksGaus(sumHist,gain,fitvals)
                if not badFit:
                    skipper_utils.fitPeaksPoisson(fitfunc,sumHist,fitvals)
                latex.DrawLatex(0.3,0.8,settingString)
            canvasNum = settingListDict[settingNames[0]].index(settingWord0) + (len(settingListDict[settingNames[1]])-1)*len(settingListDict[settingNames[0]]) + 1
            c.cd(canvasNum)
            latex.DrawLatex(0.4,0.025,settingWord0)
        for settingWord1 in settingListDict[settingNames[1]]:
            canvasNum = settingListDict[settingNames[1]].index(settingWord1)*len(settingListDict[settingNames[0]]) + 1
            c.cd(canvasNum)
            latex.SetTextAngle(90)
            latex.DrawLatex(0.05,0.5,settingWord1)
            latex.SetTextAngle(0)
        c.cd()
        latex.DrawLatex(0.4,0.95,"OHDU{0},{1}".format(ohdus[iHdu],htypes[iHist]))
        c.Print(outfilename+".pdf");

print settingListDict
c.Print(outfilename+".pdf]");
outfile.Write()
outfile.Close()
