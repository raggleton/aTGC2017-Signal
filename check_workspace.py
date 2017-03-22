from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import os

gSystem.Load("%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so"%os.environ["CMSSW_BASE"])

gSystem.Load("PDFs/PdfDiagonalizer_cc.so")
gSystem.Load("PDFs/Util_cxx.so")
gSystem.Load("PDFs/hyperg_2F1_c.so")
gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")
from ROOT import draw_error_band
from ROOT import RooErfExpPdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf


ch_nums = {"sig_el":"ch1","sb_lo_el":"ch3","sb_hi_el":"ch5","sig_mu":"ch2","sb_lo_mu":"ch4","sb_hi_mu":"ch6"}


def plot(w,spectrum,ch,region):
    
    if spectrum == "mj":
        rrv_x       = w.var("mj_%s"%region)
    elif spectrum == "mlvj":
        rrv_x       = w.var("rrv_mass_lvj")

    p           = rrv_x.frame(RooFit.Bins(rrv_x.getBins()),RooFit.Title(region+"_"+ch))
    if ch=='el':
        p.GetYaxis().SetRangeUser(1e-4,100)
    elif ch=='mu':
        p.GetYaxis().SetRangeUser(1e-4,130)

    ch_num      = ch_nums[region+'_'+ch]
    bkgs        = ['WJets','TTbar','STop']
    data        = w.data("data_obs")

    bkg_pdfs    = {}
    bkg_norms    = {}
    for bkg in bkgs:
        bkg_pdfs[bkg]       = RooExtendPdf(bkg,bkg,w.pdf("shapeBkg_%s_%s"%(bkg,ch_num)),w.function("shapeBkg_%s_%s__norm"%(bkg,ch_num)))
        bkg_norms[bkg]      = w.function("shapeBkg_%s_%s__norm"%(bkg,ch_num)).getVal()
    bkg_pdfs["WWWZ"]    = RooExtendPdf("WWWZ","WWWZ",w.pdf("shapeSig_ATGCPdf_WWWZ_%s_%s_%s"%(region,ch,ch_num)),w.function("shapeSig_ATGCPdf_WWWZ_%s_%s_%s__norm"%(region,ch,ch_num)))
    bkg_norms["WWWZ"]   = w.function("shapeSig_ATGCPdf_WWWZ_%s_%s_%s__norm"%(region,ch,ch_num)).getVal()


    model   = RooAddPdf("model","model",RooArgList(bkg_pdfs["WWWZ"],bkg_pdfs["TTbar"],bkg_pdfs["STop"],bkg_pdfs["WJets"]))
    model_norm  = float(bkg_norms["WJets"]+bkg_norms["STop"]+bkg_norms["TTbar"]+bkg_norms["WWWZ"])

    if spectrum == "mj":
        model.plotOn(p,RooFit.Name("WWWZ"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(kRed),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("TTbar"),RooFit.Components("WJets,STop,TTbar"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(kYellow),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("STop"),RooFit.Components("WJets,STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(kBlue),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("WJets"),RooFit.Components("WJets"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(kGreen),RooFit.DrawOption("F"))

        model.plotOn(p,RooFit.Name("WWWZ"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("TTbar"),RooFit.Components("WJets,STop,TTbar"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("STop"),RooFit.Components("WJets,STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(2))
        model.plotOn(p,RooFit.Name("WJets"),RooFit.Components("WJets"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(2))
    elif spectrum == "mlvj":
        model.plotOn(p,RooFit.Name("WWWZ"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(kRed),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("WJets"),RooFit.Components("STop,WJets,TTbar"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(kGreen),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("TTbar"),RooFit.Components("STop,TTbar"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(kOrange),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("STop"),RooFit.Components("STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(kBlue),RooFit.DrawOption("F"))

        model.plotOn(p,RooFit.Name("STop"),RooFit.Components("STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("TTbar"),RooFit.Components("STop,TTbar"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("WJets"),RooFit.Components("STop,WJets,TTbar"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("WWWZ"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))

    data.plotOn(p,RooFit.Cut("CMS_channel==CMS_channel::%s"%ch_num))

    return p


def plot_mj(w,ch="el"):
    
    canvas = TCanvas(ch,ch,1)

    pad_sblo    = TPad("sblo","sblo",0,0,5/22.,1)
    pad_sblo.SetRightMargin(0)
    pad_sblo.cd()
    p_sblo    = plot(w,'mj',ch,'sb_lo')
    p_sblo.GetXaxis().SetTitleSize(0)
    p_sblo.Draw()
    p_sblo.Print()

    pad_sig    = TPad("sig","sig",5/22.,0,13/22.,1)
    pad_sig.SetRightMargin(0)
    pad_sig.SetLeftMargin(0)
    pad_sig.cd()
    p_sig    = plot(w,'mj',ch,'sig')
    p_sig.GetXaxis().SetTitleSize(0)
    p_sig.Draw()

    pad_sbhi    = TPad("sbhi","sbhi",13/22.,0,1,1)
    pad_sbhi.SetLeftMargin(0)
    pad_sbhi.cd()
    p_sbhi    = plot(w,'mj',ch,'sb_hi')
    p_sbhi.Draw()

    canvas.cd()
    pad_sblo.Draw()
    pad_sig.Draw()
    pad_sbhi.Draw()

    canvas.Draw()
    canvas.Update()
    raw_input(1)

    pad_sblo.Delete()
    pad_sig.Delete()
    pad_sbhi.Delete()

def plot_mlvj(w,ch='el',reg='sig'):

    canvas = TCanvas(ch,ch,1)
    p=plot(w,'mlvj',ch,reg)
    p.GetYaxis().SetRangeUser(1e-2,2e2)
    canvas.SetLogy()
    p.Draw()
    canvas.Draw()
    canvas.Update()
    raw_input(ch)

def check_mj_4consistency(w,ch):
    bkgs    = ['WJets','TTbar','STop']
    for bkg in bkgs:
        paras_sets = []
        for region in ['sb_lo','sig','sb_hi']:
            ch_num      = ch_nums[region+'_'+ch]
            pdf     = w.pdf("shapeBkg_%s_%s"%(bkg,ch_num))
            paras   = pdf.getParameters(w.data("data_obs"))
            paras_sets.append(paras)
            
        pIter1  = paras_sets[0].createIterator()
        pIter2  = paras_sets[1].createIterator()
        pIter3  = paras_sets[2].createIterator()
        pnext1  = pIter1.Next()
        pnext2  = pIter2.Next()
        pnext3  = pIter3.Next()
        fails = []
        while pnext1:
            if pnext1!=pnext2 or pnext1!=pnext3 or pnext2!=pnext3:
                fails.extend([pnext1,pnext2,pnext3])
            pnext1 = pIter1.Next()
            pnext2 = pIter2.Next()
            pnext3 = pIter3.Next()
        print "Parameters that are different in the " + bkg + "-Pdf:"
        for i in fails:
            i.Print()
        raw_input('---')



fileIn     = TFile.Open("workspace_simfit.root")
w        = fileIn.Get("w")
fileIn.Close()
#plot_mj("el")
#plot_mj("mu")
#plot_mlvj(w,'el','sb_hi')
check_mj_4consistency(w,'mu')


