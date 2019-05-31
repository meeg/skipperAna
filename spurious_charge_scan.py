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
    ints = [int(c) for c in re.split(r'(\d+)',filename) if c.isdigit()]
    return ints[-2]

def formatVoltage(val):
    if (val<0):
        val *= -1
    if (round(val)==val):
        val = int(val)
    return str(val).replace(".","p")

outfilename = remainder[0]
infiles = remainder[1:]
infiles.sort(key=decodeRunnum)

#print(EXPOSURE,READOUT)
#if (len(sys.argv)<3):
#    print("not enough input args")
#    sys.exit(1)



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

outfile = TFile(outfilename+".root","RECREATE")

latex = TLatex()
latex.SetNDC(True)

htypes = ["allpix","prescan","overscan","active"]
regioncuts = ["x>0 && y>0","x>0 && x<8 && y>0","x>369 && y>0","x>=8 && x<=369 && y>0 && y<=624"]

numRuns = 0
settingString=""
fileHists=[]
runHists=[]
dataChain = TChain("skPixTree")
for iFile in range(0,len(infiles)+1):
    makeRunHists = False
    prevSettingString = settingString
    if (iFile!=len(infiles)):
        print(infiles[iFile])
        runnum = decodeRunnum(infiles[iFile])
        thefile = TFile(infiles[iFile])
        header = thefile.Get("headerTree_0")
        data = thefile.Get("skPixTree")
        outfile.cd()

        ogl = getHeaderValue(header,"OGAL")
        swl = getHeaderValue(header,"SWAL")
        print("OGL={0}, SWL={1}".format(ogl,swl))
        settingString = "OGL{0}SWL{1}".format(formatVoltage(ogl),formatVoltage(swl))

        if (prevSettingString!="" and settingString != prevSettingString):
            print("new settings: "+settingString)
            makeRunHists = True

    if iFile==len(infiles) or makeRunHists:
        numRuns+=1
        runHistString = prevSettingString
        histsThisRun = []
        runHists.append(histsThisRun)
        for hnum in range(0,4):
            hname = "h"+htypes[hnum]+runHistString
            dataChain.Draw("pix>>"+hname+pixbinning,regioncuts[hnum],"colz")
            h = gDirectory.Get(hname)
            histsThisRun.append(h)
        dataChain.Reset()

    if iFile==len(infiles):
        break

    dataChain.Add(infiles[iFile])

    histsThisFile = []
    fileHists.append(histsThisFile)
    fileHistString = settingString+"FILE{0}".format(runnum)

    c.Clear()
    c.Divide(1,4)
    latex.DrawLatex(0.1,0.9,"OGL={0}, SWL={1}".format(ogl,swl))

    aduZero = -1000
    aduPerElectron = -1000
    noise = -1000

    for hnum in range(0,4):
        c.cd(hnum+1)
        gPad.SetLogy(1)
        hname = "h"+htypes[hnum]+fileHistString
        data.Draw("pix>>"+hname+pixbinning,regioncuts[hnum],"colz")
        h = gDirectory.Get(hname)
        histsThisFile.append(h)
        if hnum==0:
            s0 = h.Fit("gaus","QS","",-0.4*gain,0.4*gain)
            s1 = h.Fit("gaus","QS","",0.6*gain,1.4*gain)
            badFit = (int(s0)!=0 or int(s1)!=0)
            if not badFit:
                aduZero = s0.Parameter(1)
                aduPerElectron = s1.Parameter(1) - s0.Parameter(1)
                noise = s0.Parameter(2)/aduPerElectron
                badFit = (aduPerElectron<100 or aduPerElectron>1000 or noise<0.01 or noise>1)
            if not badFit:
                s0 = h.Fit("gaus","QS","",aduZero-0.4*aduPerElectron,aduZero+0.4*aduPerElectron)
                s1 = h.Fit("gaus","QS","",aduZero+0.7*aduPerElectron,aduZero+1.4*aduPerElectron)
                badFit = (int(s0)!=0 or int(s1)!=0)
            if not badFit:
                aduZero = s0.Parameter(1)
                aduPerElectron = s1.Parameter(1) - s0.Parameter(1)
                noise = s0.Parameter(2)/aduPerElectron
                print("zero={0}, gain={1}, noise={2}".format(aduZero,aduPerElectron,noise))
                badFit = (aduPerElectron<100 or aduPerElectron>1000 or noise<0.01 or noise>1)

        if not badFit:
            fitfunc.SetParameters(h.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.01)
            h.Fit(fitfunc,"QS","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
            s = h.Fit(fitfunc,"QSL","",aduZero-0.5*aduPerElectron,aduZero+2.5*aduPerElectron)
            latex.DrawLatex(0.7,0.4,"mu={0} \pm {1}".format(s.Parameter(4),s.Error(4)))

    c.cd()
    c.Print(outfilename+".pdf");

gStyle.SetOptStat(0)
c.Clear()
c.Divide(6,4)
for iHist in range(0,4):
    for iRun in range(0,numRuns):
        h = runHists[iRun][iHist]
        settingString = h.GetName()[1+len(htypes[iHist]):]
        swlstep = iRun/4
        oglstep = iRun%4
        c.cd(oglstep*6 + swlstep + 1)
        gPad.SetLogy(1)
        h.Draw()
        latex.DrawLatex(0.3,0.8,settingString)
    c.cd()
    c.Print(outfilename+".pdf");

c.Print(outfilename+".pdf]");
outfile.Write()
outfile.Close()
