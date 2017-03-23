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


def plot(w,fitres,normset,spectrum,ch,region):
    
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

    fitparas    = fitres.floatParsFinal()
    for i in range(fitparas.getSize()):
        w.var(fitparas.at(i).GetName()).setVal(fitparas.at(i).getVal())

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

    data_histo   = data.binnedClone("data","data").createHistogram("data",rrv_x,RooFit.Cut("CMS_channel==CMS_channel::%s"%ch_num))
    data_plot    = RooHist(data_histo,100)
    data_plot.SetMarkerStyle(20)
    data_plot.SetMarkerSize(1)
    alpha        = 1-0.6827
    for iPoint in range(data_plot.GetN()):
        N = data_plot.GetY()[iPoint]
        if N==0 :
            L = 0
        else:
            L = Math.gamma_quantile(alpha/2.,N,1.)
        U = Math.gamma_quantile_c(alpha/2,N+1,1.)
        data_plot.SetPointEYlow(iPoint, N-L)
        data_plot.SetPointEYhigh(iPoint, U-N)
        data_plot.SetPointEXlow(iPoint,0)
        data_plot.SetPointEXhigh(iPoint,0)
    data_plot.SetName('data')
    p.addPlotable(data_plot,"PE")

    #data.plotOn(p,RooFit.Cut("CMS_channel==CMS_channel::%s"%ch_num),RooFit.Name("data"))

    return p


def plot_mj(w,ch="el"):
    
    canvas = TCanvas(ch,ch,1)

    pad_sblo    = TPad("sblo","sblo",0,0,5/22.,1)
    pad_sblo.SetRightMargin(0)
    pad_sblo.cd()
    p_sblo    = plot(w,fitres,normset,'mj',ch,'sb_lo')
    p_sblo.GetXaxis().SetTitleSize(0)
    p_sblo.Draw()
    p_sblo.Print()

    pad_sig    = TPad("sig","sig",5/22.,0,13/22.,1)
    pad_sig.SetRightMargin(0)
    pad_sig.SetLeftMargin(0)
    pad_sig.cd()
    p_sig    = plot(w,fitres,normset,'mj',ch,'sig')
    p_sig.GetXaxis().SetTitleSize(0)
    p_sig.Draw()

    pad_sbhi    = TPad("sbhi","sbhi",13/22.,0,1,1)
    pad_sbhi.SetLeftMargin(0)
    pad_sbhi.cd()
    p_sbhi    = plot(w,fitres,normset,'mj',ch,'sb_hi')
    p_sbhi.Draw()

    canvas.cd()
    pad_sblo.Draw()
    pad_sig.Draw()
    pad_sbhi.Draw()

    canvas.Draw()
    canvas.Update()

    raw_input(',.,.,.')

    pad_sblo.Delete()
    pad_sig.Delete()
    pad_sbhi.Delete()

def plot_mlvj(w,ch='el',reg='sig'):

    canvas = TCanvas(ch,ch,1)
    canvas.SetLogy()
    pad1        = TPad('pad1','pad1',0.,0.25,1.,1.)
    pad_pull    = TPad('pad_pull','pad_pull',0.,0.02,1.,0.25)
    p=plot(w,fitres,normset,'mlvj',ch,reg)
    p.GetYaxis().SetRangeUser(1e-2,2e2)

    canvas.cd()
    pad1.Draw()
    pad_pull.Draw()

    pad1.cd()
    pad1.SetLogy()
    pad1.SetTicky()
    p.Draw()
    p.Print()

    pad_pull.cd()
    pullhist = p.pullHist('data','WWWZ')
    ratio_style = TH1D('ratio_style','ratio_style',26,900,3500)
    ratio_style.SetMarkerStyle(21)
    ratio_style.SetMaximum(3)
    ratio_style.SetMinimum(-3)
    ratio_style.GetYaxis().SetNdivisions(7)
    ratio_style.GetYaxis().SetTitle('#frac{MC-Fit}{error}')
    ratio_style.GetYaxis().SetLabelSize(0.125)
    ratio_style.GetYaxis().SetTitleSize(0.2)
    ratio_style.GetYaxis().SetTitleOffset(0.2)
    ratio_style.Draw("")
    pullhist.SetLineColor(kBlue)
    pullhist.Draw("SAME E1")

    canvas.Update()
    raw_input(ch)

    pad1.Delete()
    pad_pull.Delete()




fileIn      = TFile.Open("workspace_simfit.root")
w           = fileIn.Get("w")
fileIn.Close()
fileIn      = TFile.Open("mlfit4plot.root")
fitres      = fileIn.Get("fit_s")
normset     = fileIn.Get("norm_fit_s")
fileIn.Close()

fitparas    = fitres.floatParsFinal()
string = 'name: pre-fit / post-fit \n'
for i in range(fitparas.getSize()):
    string += fitparas.at(i).GetName() + ': ' + str(w.var(fitparas.at(i).GetName()).getVal()) + ' / ' + str(fitparas.at(i).getVal()) + '\n'


#plot_mj(w,"mu")
plot_mlvj(w,'mu','sig')

print string



