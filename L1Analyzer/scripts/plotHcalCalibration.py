#!/usr/bin/env python
import os
import sys
import time
import logging
import argparse

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import tdrstyle as tdrstyle
tdrstyle.setTDRStyle()

logging.basicConfig(level=logging.INFO, stream=sys.stderr)

try:
    #try to load the function
    ROOT.langaufun
except AttributeError:
    #try to compile the function
    pkgdir = os.path.dirname(__file__)
    if len(pkgdir) == 0:
        pkgdir = "."
    path = os.sep.join((pkgdir, "langaus.C"))
    path = os.path.abspath(path)
    if not os.path.exists(path):
        raise Exception("ERROR: file does not exist ", path)
    ROOT.gROOT.ProcessLine(".L " + path +"+")


def calibrate(args):
    '''Calibrate the response for HCAL in Stage2Layer1'''
    ptBins = [6,9,12,15,20,25,30,35,40,45,55,70,90]

    # setup histogram
    canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

    hists = []
    hists_ptb = []
    outfile=ROOT.TFile(args.inputFile)
    for i in range(0,len(ptBins)):
        hname_ptb = 'hist_ptb_{0}'.format(i)
        hists_ptb.append(ROOT.TH1F(hname_ptb,'',28,0,28))
        #hists_ptb.append(outfile.Get(hname_ptb))
        hists.append([])
        #for j in range(0,56):
        for j in range(0,28):
            hname = 'hist_{0}_{1}'.format(i, j)
            hists[i].append(outfile.Get(hname))

    
    # calculate calibration
    colors = [
        ROOT.kBlack,
        ROOT.kGreen+3,
        ROOT.kBlue+2,
        ROOT.kOrange,
        ROOT.kGreen-3,
        ROOT.kPink-3,
        ROOT.kBlue-6,
        ROOT.kYellow+3,
        ROOT.kAzure+6,
        ROOT.kRed-7,
        ROOT.kMagenta+2,
        ROOT.kMagenta+3,
        ROOT.kMagenta+4,
    ]

    with open("h.txt","w") as text_file:
        for ptb in range(0,len(ptBins)):
            text_file.write("\n")
            #for eta in range(0,56):
            for eta in range(0,28):
                if not hists[ptb][eta]:
                    print 'No hist',ptb,eta
                    val = 1.
                    err = 0.
                else:
                    val = hists[ptb][eta].GetMean()
                    err = hists[ptb][eta].GetMeanError()
                    #hists[i][j].Fit('gaus')
                    #hists[i][j].Fit('landau')
                    tf1 = ROOT.TF1("landaugausfunction", ROOT.langaufun, 0, 5, 4)
                    rms = hists[ptb][eta].GetRMS()
                    peakpos = hists[ptb][eta].GetXaxis().GetBinCenter(hists[ptb][eta].GetMaximumBin())
                    startwidth = rms / 5.0
                    startmpv = peakpos
                    startnorm = hists[ptb][eta].Integral()
                    startsigma = rms / 10.0
                    tf1.SetParNames("LandauWidth","LandauMPV","Normalisation","GaussianSigma")
                    tf1.SetParameters(startwidth, startmpv, startnorm, startsigma)
                    print 'start parameters', ptb, eta
                    print startwidth, startmpv, startnorm, startsigma
                    hists[ptb][eta].Fit(tf1, "0L", "", 0, 5)
                    hists[ptb][eta].Draw()
                    tf1.Draw('same')
                    val = float(tf1.GetParameter(1))
                    err = float(tf1.GetParError(1))
                    canvas.SaveAs('hist_{0}_{1}.png'.format(ptb,eta))
                    if val!=val:
                        val = 1.
                        err = 0.
                print 'fit', ptb, eta, val, err
                text_file.write("%f, " % val)
                etab = eta+1
                hists_ptb[ptb].SetBinContent(etab,val)
                hists_ptb[ptb].SetBinError(etab,err)
        legend = ROOT.TLegend(0.25,0.8,0.85,0.9,'','brNDC')
        legend.SetNColumns(3)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetBorderSize(0)
        for ptb in range(0,len(ptBins)):
            hists_ptb[ptb].SetMarkerColor(colors[ptb])
            hists_ptb[ptb].SetLineColor(colors[ptb])
            hists_ptb[ptb].SetMarkerStyle(23)
            if ptb==0:
                hists_ptb[ptb].GetXaxis().SetTitle("TPG iEta")
                hists_ptb[ptb].GetYaxis().SetTitle("RECO p_{T} / TPG p_{T}")
                hists_ptb[ptb].GetYaxis().SetRangeUser(0.,3.1)
                #hists_ptb[ptb].GetXaxis().SetRangeUser(-28,28)
                hists_ptb[ptb].GetXaxis().SetRangeUser(0,28)
                hists_ptb[ptb].SetTitle('')
                hists_ptb[ptb].Draw('pE1')
                legend.AddEntry(hists_ptb[ptb], '3 < p_{{T}} < {0} GeV'.format(ptBins[ptb]))
            else:
                hists_ptb[ptb].Draw('pE1same')
                legend.AddEntry(hists_ptb[ptb], '{0} < p_{{T}} < {1} GeV'.format(ptBins[ptb-1],ptBins[ptb]))
        legend.Draw('same')
        savename = 'SF_HCAL.png'
        canvas.SaveAs(savename) 

    #with open('h.txt','r') as f, open('hplus.txt','w') as fp, open('hminus.txt','w') as fm:
    #     for row in f.readlines():
    #         vals = [x.strip() for x in row.split(',')]
    #         if len(vals)<56: continue
    #         fp.write(', '.join(vals[28:56])+',\n')
    #         fm.write(', '.join(reversed(vals[0:28]))+',\n')

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Validate Ntuples')

    parser.add_argument('inputFile', type=str, help='Input flat file')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    calibrate(args)

if __name__ == "__main__":
    status = main()
    sys.exit(status)
