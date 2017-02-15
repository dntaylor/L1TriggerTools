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



def calibrate(args):
    '''Calibrate the response for ECAL in Stage2Layer1'''
    # setup histogram
    canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)
    hists = []
    hists_ptb = []
    for i in range(0,9):
        hname_ptb = 'hist_ptb_{0}'.format(i)
        hists_ptb.append(ROOT.TH1F(hname_ptb,'',56,-28,28))
        hists.append([])
        for j in range(0,56):
            hname = 'hist_{0}_{1}'.format(i, j)
            hists[i].append( ROOT.TH1F(hname,'',100,0,5) )

    # fill
    collectionName = 'electron'
    for fname in args.inputFiles:
        tfile = ROOT.TFile.Open(fname)
        tree = tfile.Get(args.treeName)
        for row in tree:
            n = len(getattr(row,'{0}_pt'.format(collectionName)))
            for i in range(n):
                pt = getattr(row,'{0}_pt'.format(collectionName))[i]
                tp_et = getattr(row,'{0}_matchedEt'.format(collectionName))[i]
                sumEt = getattr(row,'{0}_matchedSum3x3'.format(collectionName))[i]
                ecalSumEt = getattr(row,'{0}_matchedEcalSum3x3'.format(collectionName))[i]
                if pt>5 and ecalSumEt>5 and tp_et/sumEt>0.25:
                    sf = pt/ecalSumEt
                    ptBin = int(getattr(row,'{0}_matchedPtBin'.format(collectionName))[i])
                    ieta = getattr(row,'{0}_matchedIEta'.format(collectionName))[i]
                    etaBin = int(ieta+28)
                    if ieta>0: etaBin-=1
                    hists[ptBin][etaBin].Fill(sf)

    # calculate calibration
    colors = [ROOT.kBlue+2, ROOT.kOrange, ROOT.kGreen-3, ROOT.kPink-3, ROOT.kBlue-6, ROOT.kYellow+3, ROOT.kAzure-6, ROOT.kRed-7, ROOT.kMagenta+2]
    with open("eg.txt","w") as text_file:
        legend = ROOT.TLegend(0.2,0.6,0.45,0.9,'','brNDC')
        legend.SetFillColor(ROOT.kWhite)
        legend.SetBorderSize(0)
        for ptb in range(0,9):
            text_file.write("\n")
            for eta in range(0,56):
                Mean = hists[ptb][eta].GetMean()
                MeanError = hists[ptb][eta].GetMeanError()
                text_file.write("%f, " % Mean)
                etab = eta+1
                hists_ptb[ptb].SetBinContent(etab,Mean)
                hists_ptb[ptb].SetBinError(etab,MeanError)
            hists_ptb[ptb].SetMarkerColor(colors[ptb])
            hists_ptb[ptb].SetLineColor(colors[ptb])
            hists_ptb[ptb].SetMarkerStyle(23)
            if ptb==0:
                hists_ptb[ptb].GetXaxis().SetTitle("TPG iEta")
                hists_ptb[ptb].GetYaxis().SetTitle("RECO p_{T} / TPG p_{T}")
                hists_ptb[ptb].GetYaxis().SetRangeUser(0.7,2.1)
                hists_ptb[ptb].GetXaxis().SetRangeUser(-28,28)
                hists_ptb[ptb].SetTitle('')
                hists_ptb[ptb].Draw('pE1')
                legend.AddEntry(hists_ptb[ptb], 'p_{{T}} < {0} GeV'.format(5*ptb+10))
            elif ptb==8:
                hists_ptb[ptb].Draw('pE1same')
                legend.AddEntry(hists_ptb[ptb], 'p_{{T}} > {0} GeV'.format(5*ptb+5))
            else:
                hists_ptb[ptb].Draw('pE1same')
                legend.AddEntry(hists_ptb[ptb], '{0} < p_{{T}} < {1} GeV'.format(5*ptb+5, 5*ptb+10))
        legend.Draw('same')
        savename = 'SF_ECAL.png'
        canvas.SaveAs(savename) 

    outfile=ROOT.TFile("outfile.root","RECREATE")
    outfile.cd()
    
    for ptb in range(0,9):
        hists_ptb[ptb].Write()
        for eta in range(0,28):
            hists[ptb][eta].Write() 

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Validate Ntuples')

    parser.add_argument('inputFiles', nargs='*', help='List of files with calibration tree')
    parser.add_argument('--treeName', type=str, default='l1Analyzer/L1TTree', help='Name of tree')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    calibrate(args)

if __name__ == "__main__":
    status = main()
    sys.exit(status)
