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
    '''Calibrate the response for HCAL in Stage2Layer1'''
    ptBins = [6,9,12,15,20,25,30,35,40,45,55,70,90]

    # setup histogram
    canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)
    hists = []
    hists_ptb = []
    for i in range(0,len(ptBins)):
        hname_ptb = 'hist_ptb_{0}'.format(i)
        #hists_ptb.append(ROOT.TH1D(hname_ptb,'',56,-28,28))
        hists_ptb.append(ROOT.TH1D(hname_ptb,'',28,0,28))
        hists.append([])
        #for j in range(0,56):
        for j in range(0,28):
            hname = 'hist_{0}_{1}'.format(i, j)
            hists[i].append( ROOT.TH1D(hname,'',100,0,5) )


    # fill
    collectionName = 'genPion'
    for f,fname in enumerate(args.inputFiles):
        logging.info('Processing file {0} of {1}'.format(f+1, len(args.inputFiles)))
        tfile = ROOT.TFile.Open(fname)
        tree = tfile.Get(args.treeName)
        for row in tree:
            n = len(getattr(row,'{0}_pt'.format(collectionName)))
            for i in range(n):
                iphi = getattr(row,'{0}_matchedIPhi'.format(collectionName))[i]
                ieta = getattr(row,'{0}_matchedIEta'.format(collectionName))[i]
                if args.hep17:
                    if iphi<63 or iphi>66 or ieta<0: continue
                if args.nothep17:
                    if iphi>=60 and iphi<=69 and ieta>0: continue
                pt = getattr(row,'{0}_pt'.format(collectionName))[i]
                tp_et = getattr(row,'{0}_matchedEt'.format(collectionName))[i]
                if args.closure:
                    tp_et = getattr(row,'{0}_matchedCorrectedEt'.format(collectionName))[i]
                ecalSumEt = getattr(row,'{0}_matchedEcalCorrectedSum5x5'.format(collectionName))[i]
                hcalSumEt = getattr(row,'{0}_matchedHcalSum5x5'.format(collectionName))[i]
                if args.closure:
                    hcalSumEt = getattr(row,'{0}_matchedHcalCorrectedSum5x5'.format(collectionName))[i]
                #ecalSumEt3x3 = getattr(row,'{0}_matchedEcalCorrectedSum3x3'.format(collectionName))[i]
                #hcalSumEt3x3 = getattr(row,'{0}_matchedHcalSum3x3'.format(collectionName))[i]
                sumEt = ecalSumEt+hcalSumEt
                #sumEt3x3 = ecalSumEt3x3+hcalSumEt3x3
                #if pt>5 and hcalSumEt>5 and tp_et/sumEt>0.25:
                #if pt>5 and hcalSumEt>5 and tp_et/sumEt>0.95:
                #if abs(ieta)>20 and abs(ieta)<29 and pt>5 and hcalSumEt>5:
                #    print ieta, iphi, pt, tp_et, ecalSumEt, hcalSumEt, sumEt, hcalSumEt3x3
                if abs(ieta)>20: tp_et *= 2
                #if pt>5 and hcalSumEt>3 and pt<128 and tp_et/hcalSumEt>0.75:
                #if pt>0 and hcalSumEt>3 and tp_et/hcalSumEt>0.95:
                if pt>0 and hcalSumEt>3 and tp_et/hcalSumEt>0.2 and tp_et<90:
                    sf = pt/sumEt
                    #sf = (pt-ecalSumEt)/hcalSumEt
                    ptBin = 0
                    for i,high in enumerate(ptBins):
                        if tp_et<high: break
                        ptBin = i+1
                    #etaBin = int(ieta+28)
                    #if ieta>0: etaBin-=1
                    #if etaBin>56 or etaBin<0: continue # HF
                    etaBin = int(abs(ieta)-1)
                    if etaBin>28:continue
                    hists[ptBin][etaBin].Fill(sf)

    outfile=ROOT.TFile(args.outputFile,"RECREATE")
    outfile.cd()
    
    for ptb in range(0,len(ptBins)):
        hists_ptb[ptb].Write()
        #for eta in range(0,56):
        for eta in range(0,28):
            hists[ptb][eta].Write() 

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Validate Ntuples')

    parser.add_argument('inputFiles', nargs='*', help='List of files with calibration tree')
    parser.add_argument('--treeName', type=str, default='l1Analyzer/L1TTree', help='Name of tree')
    parser.add_argument('--outputFile', type=str, default='outfile_hcal.root', help='Output file name')
    parser.add_argument('--closure', action='store_true', help='Closure test')
    parser.add_argument('--hep17', action='store_true', help='Select HEP17')
    parser.add_argument('--nothep17', action='store_true', help='Select not HEP17')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    calibrate(args)

if __name__ == "__main__":
    status = main()
    sys.exit(status)
