from ROOT import  *
from array import array
from optparse import OptionParser
from ConfigParser import SafeConfigParser
import math as math
import random
import os

gSystem.Load('%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so'%os.environ['CMSSW_BASE'])
from ROOT import RooErfExpPdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf
gSystem.Load('PDFs/PdfDiagonalizer_cc.so')
gSystem.Load('PDFs/Util_cxx.so')
gSystem.Load('PDFs/hyperg_2F1_c.so')
gSystem.Load('PDFs/HWWLVJRooPdfs_cxx.so')
from ROOT import draw_error_band


parser        = OptionParser()
parser.add_option('-n', '--newtrees', action='store_true', dest='newtrees', default=False, help='recreate aTGC histograms')
parser.add_option('-p', '--plots', action='store_true', dest='Make_plots', default=False, help='make plots')
parser.add_option('--savep', action='store_true', dest='savep', default=False, help='save plots')
parser.add_option('--p2', action='store_true', dest='Make_plots2', default=False, help='make parabel plot')
parser.add_option('-b', action='store_true', dest='batch', default=False, help='batch mode')
parser.add_option('-c', '--ch', dest='chan', default='elmu', help='channel, el, mu or elmu')
parser.add_option('--noatgcint', action='store_true', dest='noatgcint', default=False, help='set atgc-interference coefficients to zero')
parser.add_option('--printatgc', action='store_true', default=False, help='print atgc-interference contribution')
parser.add_option('--yieldplots', action='store_true', default=False, help='make plots of relative yields')
parser.add_option('--atgc', action='store_true', dest='atgc', default=False, help='use anomalous coupling parametrization instead of EFT')
parser.add_option("--config",dest="config",default='config_simfit',help="The name of the input configuration file.")


(options,args) = parser.parse_args()


gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
if options.batch:
        gROOT.SetBatch(True)

if not os.path.isdir('Output'):
        os.system('mkdir Output')
if not os.path.isdir('docuplots'):
        os.system('mkdir docuplots')


class Prepare_workspace_4limit:
    
        def __init__(self,ch,mlvj_lo,mlvj_hi):
            
            self.POI                    = ['cwww','ccw','cb']
            self.PAR_TITLES             = {'cwww' : '#frac{c_{WWW}}{#Lambda^{2}}', 'ccw' : '#frac{c_{W}}{#Lambda^{2}}', 'cb' : '#frac{c_{B}}{#Lambda^{2}}'}#latex titles 
            self.PAR_MAX                = {'cwww' : 12, 'ccw' : 20, 'cb' : 60}#atgc points
            self.signal_category        = 'WV'#only one signal region
            self.ch                     = ch
            
            self.mj_lo                  = 65.                    #lower bound of the signal region
            self.mj_hi                  = 105.                   #upper bound
            self.binlo                  = mlvj_lo                #lower bound on invariant mass
            self.binhi                  = mlvj_hi                #upper bound

            self.channel                = "WV_"+self.ch
            self.nbins                  = (self.binhi-self.binlo)/100
            
            self.WS                     = RooWorkspace("w")        #final workspace
            self.wtmp                   = RooWorkspace('wtmp')
            
            self.fitresults             = []
            ##nuisance parameter to change all slope parameters by certain percentage (bigger for cb in WZ-cateogry)
            self.eps                    = RooRealVar('slope_nuis','slope_nuis',1,0,2)
            self.eps.setConstant(kTRUE)
            self.eps4cbWZ               = RooFormulaVar('rel_slope_nuis4cbWZ','rel_slope_nuis4cbWZ','1+3*(@0-1)',RooArgList(self.eps))
            
            ##read workspace containing background pdfs
            fileInWs                    = TFile.Open('Input/wwlvj_%s_HPV_workspace.root'%(self.ch))
            w                           = fileInWs.Get('workspace4limit_')
            self.rrv_mass_lvj           = w.var('rrv_mass_lvj')
            self.rrv_mass_lvj.SetTitle('M_{WV}')
            self.rrv_mass_lvj.setRange(self.binlo,self.binhi)
            self.rrv_mass_j             = w.var('rrv_mass_j')

            #names of bkg contributions and regions
            self.bkgs        = ['WJets','TTbar','WW','WZ','STop']
            self.regions     = ['sig','sb_lo','sb_hi']

            #read config##FIXME adjust naming
            self.cfg = SafeConfigParser()
            self.cfg.read(options.config)


        #read trees containing aTGC WW and WZ events and fill them into histograms
        def Read_ATGCtree(self,ch='el'):

            #cut on MET has to be applied
            if self.ch=='el':
                METCUT        = 80.
            elif self.ch=='mu':
                METCUT        = 40.
            else:
                raise RuntimeError('no such channel %s'%self.ch)
            
            print '##########making histograms for aTGC working points##########'
            hists4scale        = {}
            
            for WV in ['WW','WZ']:
                #create 3 histograms for each aTGC parameter (positive, negative and positive-negative working point)
                for para in self.POI:
                    hists4scale['c_pos_%s_hist_%s'%(WV,para)] = TH1F('c_pos_%s_hist_%s'%(WV,para),'c_pos_%s_hist_%s'%(WV,para),self.nbins,self.binlo,self.binhi);
                    hists4scale['c_neg_%s_hist_%s'%(WV,para)] = TH1F('c_neg_%s_hist_%s'%(WV,para),'c_neg_%s_hist_%s'%(WV,para),self.nbins,self.binlo,self.binhi);
                    hists4scale['c_dif_%s_hist_%s'%(WV,para)] = TH1F('c_dif_%s_hist_%s'%(WV,para),'dif_%s_hist_%s'%(WV,para),self.nbins,self.binlo,self.binhi);
                #add histograms for SM and all aTGC parameters unequal to zero
                hists4scale['c_SM_%s_hist'%WV]                  = TH1F('c_SM_%s_hist'%WV,'c_SM_%s_hist'%WV,self.nbins,self.binlo,self.binhi);                
                hists4scale['c_%s_histall3'%WV]                 = TH1F('c_%s_histall3'%WV,'c_%s_histall3'%WV,self.nbins,self.binlo,self.binhi);

                print 'reading %s-aTGC_%s.root'%(WV,self.ch)
                fileInATGC        = TFile.Open('Input/%s-aTGC_%s.root'%(WV,self.ch))
                treeInATGC        = fileInATGC.Get('treeDumper/BasicTree')

                lumi_tmp         = 2300.

                for i in range(treeInATGC.GetEntries()):
                    if i%10000==0:
                            print str(i) + '/' + str(treeInATGC.GetEntries())
                    treeInATGC.GetEntry(i)
                    MWW                = treeInATGC.MWW
                    #apply cuts
                    if treeInATGC.jet_pt>200. and treeInATGC.jet_tau2tau1<0.6 and treeInATGC.W_pt>200. and treeInATGC.deltaR_LeptonWJet>math.pi/2. and treeInATGC.Mjpruned>40 and treeInATGC.Mjpruned<150 and abs(treeInATGC.deltaPhi_WJetMet)>2. and abs(treeInATGC.deltaPhi_WJetWlep)>2. and treeInATGC.nbtag==0 and treeInATGC.pfMET>METCUT and MWW>self.binlo:
                        weight_part = 1/20. * lumi_tmp * treeInATGC.totWeight
                        aTGC        = treeInATGC.aTGCWeights                #contains weights for different workingpoints
                        #all3
                        hists4scale['c_%s_histall3'%WV].Fill(MWW,aTGC[6] * weight_part)
                        #SM
                        hists4scale['c_SM_%s_hist'%WV].Fill(MWW,aTGC[7] * weight_part)
                        #cwww
                        hists4scale['c_pos_%s_hist_cwww'%WV].Fill(MWW,aTGC[0] * weight_part)
                        hists4scale['c_neg_%s_hist_cwww'%WV].Fill(MWW,aTGC[1] * weight_part)
                        #ccw
                        hists4scale['c_pos_%s_hist_ccw'%WV].Fill(MWW,aTGC[2] * weight_part)
                        hists4scale['c_neg_%s_hist_ccw'%WV].Fill(MWW,aTGC[3] * weight_part)
                        #cb
                        hists4scale['c_pos_%s_hist_cb'%WV].Fill(MWW,aTGC[4] * weight_part)
                        hists4scale['c_neg_%s_hist_cb'%WV].Fill(MWW,aTGC[5] * weight_part)
                        #ccw-SM interference
                        hists4scale['c_dif_%s_hist_ccw'%WV].Fill(MWW,aTGC[2]-aTGC[3])
                        #cb-SM interference
                        hists4scale['c_dif_%s_hist_cb'%WV].Fill(MWW,aTGC[4]-aTGC[5])
                            
            #write histograms to file
            fileOut        = TFile.Open('Output/hists4scale_%s_WV_aTGC-%s_%s.root'%(self.ch,self.binlo,self.binhi),'recreate')
            for key in hists4scale:
                hists4scale[key].Write()
            print '--------> Written to file ' + fileOut.GetName()
            fileOut.Close()

        def Make_plots(self,rrv_x,cat,fitres):
                
            can     = []
            can2    = []
            plots   = []
            plots2  = []
            pads    = []
            
            channel = self.ch+'_'+cat

            for i in range(3):
                rrv_x.setRange(self.binlo,self.binhi)
                p       = rrv_x.frame(self.binlo,self.binhi)
                p2      = rrv_x.frame(self.binlo,self.binhi)
                c       = TCanvas(self.POI[i]+'-',self.POI[i]+'-',1)
                c.cd()
                pad1        = TPad('pad1_%s'%self.POI[i],'pad1_%s'%self.POI[i],0.,0.25,1.,1.)
                pad2        = TPad('pad2_%s'%self.POI[i],'pad2_%s'%self.POI[i],0.,0.02,1.,0.25)
                c2          = TCanvas(self.POI[i]+'+',self.POI[i]+'+',1)
                c2.cd()
                pad3        = TPad('pad3_%s'%self.POI[i],'pad3_%s'%self.POI[i],0.,0.25,1.,1.)
                pad4        = TPad('pad4_%s'%self.POI[i],'pad4_%s'%self.POI[i],0.,0.02,1.,0.25)
                p2pads      = [pad1,pad2,pad3,pad4]
                can.append(c)
                can2.append(c2)
                plots.append(p)
                plots2.append(p2)
                pads.append(p2pads)

            for i in range(3):

                can[i].cd()
                pads[i][0].Draw()
                pads[i][1].Draw()
                pads[i][0].SetLeftMargin(0.1)
                pads[i][1].SetLeftMargin(0.1)
                
                norm = self.wtmp.function('normfactor_3d_%s'%channel)

                for j in range(3):
                        self.wtmp.var(self.POI[j]).setVal(0)
                self.wtmp.data('SMdatahist_%s'%cat).plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'),RooFit.Name('SMdata'))
                normvalSM        = norm.getVal() * self.wtmp.data('SMdatahist_%s'%cat).sumEntries()

                self.wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots[i],RooFit.LineColor(kBlack),RooFit.Normalization(normvalSM, RooAbsReal.NumEvent),RooFit.Name('SMmodel'))

                self.wtmp.data('neg_datahist_%s_%s'%(cat,self.POI[i])).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'),RooFit.Name('atgcdata'))
                self.wtmp.var(self.POI[i]).setVal(-self.PAR_MAX[self.POI[i]])
                normvalneg = norm.getVal() * self.wtmp.data('SMdatahist_%s'%cat).sumEntries()
                self.wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots[i],RooFit.LineColor(kBlue),RooFit.Normalization(normvalneg, RooAbsReal.NumEvent),RooFit.Name('atgcmodel'))
        
                pullhist = plots[i].pullHist('atgcdata','atgcmodel')
                
                plotmax        = 100
                if self.ch == 'el':
                    plotmin = 1e-4
                    if cat == 'WZ':
                        plotmin = 3e-4
                elif self.ch == 'mu':
                    plotmin = 1e-2
                    if cat == 'WZ':
                        plotmin = 1e-3
                        plotmax = 50
                if cat == 'WV':
                    plotmin = 1e-2
                    plotmax = 1.5e2
                plots[i].GetYaxis().SetRangeUser(plotmin,plotmax)
                pads[i][0].cd()
                pads[i][0].SetLogy()
                plots[i].SetTitle('')
                plots[i].GetYaxis().SetTitle('arbitrary units')
                #plots[i].GetXaxis().SetRangeUser(900,3500)

                plots[i].Draw()
                ndof        = (self.binhi-self.binlo)/100 - 4
                plots[i].Print()
                #print calculate_chi2(self.wtmp.data('neg_datahist4fit_%s'%self.POI[i]),rrv_x,pullhist,plots[i],ndof,1)
                
                parlatex        = ['#frac{c_{WWW}}{#Lambda^{2}}','#frac{c_{W}}{#Lambda^{2}}','#frac{c_{B}}{#Lambda^{2}}']
                leg        = TLegend(0.11,0.2,0.4,0.6)
                leg.SetFillStyle(0)
                leg.SetBorderSize(0)
                leg.AddEntry(plots[i].findObject('SMdata'),'MC '+parlatex[i]+'=0 TeV^{-2}','le')
                leg.AddEntry(plots[i].findObject('SMmodel'),'signal model '+parlatex[i]+'=0 TeV^{-2}','l')
                leg.AddEntry(plots[i].findObject('atgcdata'),'MC '+parlatex[i]+'='+str(-self.PAR_MAX[self.POI[i]])+' TeV^{-2}','le')
                leg.AddEntry(plots[i].findObject('atgcmodel'),'signal model '+parlatex[i]+'='+str(-self.PAR_MAX[self.POI[i]])+' TeV^{-2}','l')
                leg.Draw()
                leg.Print()
                
                pads[i][1].cd()
                ratio_style = TH1D('ratio_style','ratio_style',(self.binhi-self.binlo)/100,self.binlo,self.binhi)
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


                can[i].Update()
                if options.savep:
                    can[i].SaveAs('tmpplots/%s_neg_%s.pdf'%(self.POI[i],channel))
                    can[i].SaveAs('tmpplots/%s_neg_%s.png'%(self.POI[i],channel))
                

                for j in range(3):
                        self.wtmp.var(self.POI[j]).setVal(0)
                self.wtmp.data('SMdatahist_%s'%cat).plotOn(plots2[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'))
                self.wtmp.data('pos_datahist_%s_%s'%(cat,self.POI[i])).plotOn(plots2[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'))

                self.wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],RooFit.LineColor(kBlack),RooFit.Normalization(normvalSM, RooAbsReal.NumEvent))
                self.wtmp.var(self.POI[i]).setVal(self.PAR_MAX[self.POI[i]])
                normvalpos = norm.getVal() * self.wtmp.data('SMdatahist_%s'%cat).sumEntries()

                self.wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],RooFit.LineColor(kBlue),RooFit.Normalization(normvalpos, RooAbsReal.NumEvent))
                
                self.wtmp.data('pos_datahist_%s_%s'%(cat,self.POI[i])).plotOn(plots2[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
                plots2[i].GetYaxis().SetRangeUser(plotmin,plotmax)
                plots2[i].GetYaxis().SetTitle('arbitrary units')
                can2[i].cd()
                pads[i][2].Draw()
                pads[i][3].Draw()
                pads[i][2].SetLeftMargin(0.1)
                pads[i][3].SetLeftMargin(0.1)
                plots2[i].SetTitle('')
                pads[i][2].SetLogy()
                pads[i][2].cd()
                plots2[i].Draw()
                pullhist2 = plots2[i].pullHist('h_pos_datahist_%s_%s'%(cat,self.POI[i]),'aTGC_model_%s_Norm[rrv_mass_lvj]'%channel)
                pads[i][3].cd()
                ratio_style.Draw("")
                pullhist2.SetLineColor(kBlue)
                pullhist2.Draw("E1")

                can2[i].Update()
                if options.savep:
                    can2[i].SaveAs('tmpplots/%s_pos_%s.pdf'%(self.POI[i],channel))
                    can2[i].SaveAs('tmpplots/%s_pos_%s.png'%(self.POI[i],channel))
                    
            if not options.batch:
                raw_input('plots plotted')

            
        #function to import multiple items from a list into a workspace
        def Import_to_ws(self,workspace,items,recycle=0):
            for item in items:
                if recycle:
                    getattr(workspace,'import')(item,RooFit.RecycleConflictNodes())
                else:
                    getattr(workspace,'import')(item)


        def Make_signal_pdf(self,rrv_x,sample):
            
            channel        = self.ch+'_'+sample                #needed for variables that differ for WW and WZ
            
            cwww    = RooRealVar('cwww','cwww',0,-120,120);
            ccw     = RooRealVar('ccw','ccw',0,-200,200);
            cb      = RooRealVar('cb','cb',0,-600,600);
            cwww.setConstant(kTRUE);
            ccw.setConstant(kTRUE);
            cb.setConstant(kTRUE);
   
            #get SM histogram and make RooDataHist
            fileInHist      = TFile.Open('Output/hists4scale_%s_WV_aTGC-%s_%s.root'%(self.ch,self.binlo,self.binhi))
            rrv_x.setRange(self.binlo,self.binhi)
            SMdatahist      = RooDataHist('SMdatahist_%s'%sample,'SMdatahist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_SM_%s_hist'%sample))
            fileInHist.Close()

            #make SM pdf, simple exponential
            a1_4fit         = RooRealVar('a_SM_4fit_%s'%channel,'a_SM_4fit_%s'%channel,-0.1,-2,0)
            a1              = RooFormulaVar('a_SM_%s'%channel,'a_SM_%s'%channel,'@0*@1',RooArgList(a1_4fit,self.eps))
            SMPdf           = RooExponential('SMPdf_%s'%channel,'SMPdf_%s'%channel,rrv_x,a1)
            ##actual fit to determine SM shape parameter a1_4fit
            fitresSM        = SMPdf.fitTo(SMdatahist, RooFit.SumW2Error(kTRUE), RooFit.Save(kTRUE))
            self.fitresults.append(fitresSM)
            a1_4fit.setConstant(kTRUE)
            #coefficient for SM term in final signal function
            N_SM            = RooRealVar('N_SM_%s'%channel,'N_SM_%s'%channel,SMdatahist.sumEntries())
            N_SM.setConstant(kTRUE)

            self.Import_to_ws(self.wtmp,[cwww,ccw,cb,self.eps4cbWZ,SMdatahist,SMdatahist,N_SM])
            
            #define parameter ranges for error function, only needed for WZ
            if sample=='WZ':
                if self.ch=='el':
                    Erf_width_cwww      = RooRealVar('Erf_width_cwww_%s'%channel,'Erf_width_cwww_%s'%channel,1000.,500.,1500.)
                    Erf_width_ccw       = RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,1500.,1000.,2000.)
                elif self.ch=='mu':
                    Erf_width_cwww      = RooRealVar('Erf_width_cwww_%s'%channel,'Erf_width_cwww_%s'%channel,1000.,500.,7500.)
                    Erf_width_ccw       = RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,1500.,500.,2000.)
                Erf_offset_cwww     = RooRealVar('Erf_offset_cwww_%s'%channel,'Erf_offset_cwww_%s'%channel,1000.,500.,1500.)
                Erf_offset_ccw      = RooRealVar('Erf_offset_ccw_%s'%channel,'Erf_offset_ccw_%s'%channel,1500.,500.,2500.)
                Erf_offset_cwww.setConstant(kTRUE)
                Erf_width_cwww.setConstant(kTRUE)
                Erf_offset_ccw.setConstant(kTRUE)
                Erf_width_ccw.setConstant(kTRUE)
                self.Import_to_ws(self.wtmp,[Erf_width_cwww,Erf_offset_cwww,Erf_width_ccw,Erf_offset_ccw])
                
            for i in range(len(self.POI)):
                s_name          = self.POI[i] + '_' + channel #added to parameter names
                fileInHist      = TFile.Open('Output/hists4scale_%s_WV_aTGC-%s_%s.root'%(self.ch,self.binlo,self.binhi))
                rrv_x.setRange(self.binlo,self.binhi)                
                pos_datahist    = RooDataHist('pos_datahist_%s_%s'%(sample,self.POI[i]),'pos_datahist_%s_%s'%(sample,self.POI[i]),RooArgList(rrv_x),fileInHist.Get('c_pos_%s_hist_%s'%(sample,self.POI[i])))
                neg_datahist    = RooDataHist('neg_datahist_%s_%s'%(sample,self.POI[i]),'neg_datahist_%s_%s'%(sample,self.POI[i]),RooArgList(rrv_x),fileInHist.Get('c_neg_%s_hist_%s'%(sample,self.POI[i])))
                dif_datahist    = RooDataHist('dif_datahist_%s_%s'%(sample,self.POI[i]),'dif_datahist_%s_%s'%(sample,self.POI[i]),RooArgList(rrv_x),fileInHist.Get('c_dif_%s_hist_%s'%(sample,self.POI[i])))
                fileInHist.Close()
                
                #import datasets to wtmp and final workspace WS
                self.Import_to_ws(self.wtmp,[pos_datahist,neg_datahist,dif_datahist])
                self.Import_to_ws(self.WS,[pos_datahist,neg_datahist,dif_datahist])
                
                #get scaling parabel from yields
                hist4scale = TH1F('hist4scale_%s'%self.POI[i],'hist4scale_%s'%self.POI[i],3,-1.5*self.PAR_MAX[self.POI[i]],1.5*self.PAR_MAX[self.POI[i]])
                hist4scale.SetBinContent(1,neg_datahist.sumEntries()/SMdatahist.sumEntries())
                hist4scale.SetBinContent(2,1)
                hist4scale.SetBinContent(3,pos_datahist.sumEntries()/SMdatahist.sumEntries())
                #fit parabel
                hist4scale.Fit('pol2','0')
                fitfunc     = hist4scale.GetFunction('pol2')
                par1        = RooRealVar('par1_%s'%s_name,'par1_%s'%s_name,fitfunc.GetParameter(1));
                par1.setConstant(kTRUE);
                par2        = RooRealVar('par2_%s'%s_name,'par2_%s'%s_name,fitfunc.GetParameter(2));
                par2.setConstant(kTRUE);

                N_pos_tmp   = pos_datahist.sumEntries()
                N_neg_tmp   = neg_datahist.sumEntries()
                N_quad      = RooRealVar('N_quad_%s'%s_name,'N_quad_%s'%s_name, ((N_pos_tmp+N_neg_tmp)/2)-N_SM.getVal() )
                
                #scaleshape is the relative change to SM
                scaleshape  = RooFormulaVar('scaleshape_%s'%s_name,'scaleshape_%s'%s_name, '(@0*@2+@1*@2**2)', RooArgList(par1,par2,self.wtmp.var(self.POI[i])))
                #FIXME only very few atgc events for cb in WZ sample, fit doesn't work yet -> different parametrization, or leave out completely
                if sample=='WZ' and self.POI[i]=='cb':
                    N_lin       = RooRealVar('N_lin_%s'%s_name,'N_lin_%s'%s_name, 0)
                    a2_4fit     = RooRealVar('a_quad_4fit_%s'%s_name,'a_quad_4fit_%s'%s_name,-0.1,-2,0.)
                    a2          = RooFormulaVar('a_quad_nuis_%s'%s_name,'a_quad_nuis_%s'%s_name,'@0*@1',RooArgList(a2_4fit,self.eps4cbWZ))
                    a3_4fit     = RooRealVar('a_lin_4fit_%s'%s_name,'a_lin_4fit_%s'%s_name,-0.0001,-0.1,0.)
                    a3          = RooFormulaVar('a_lin_nuis_%s'%s_name,'a_lin_nuis_%s'%s_name,'@0*@1',RooArgList(a3_4fit,self.eps4cbWZ))
                    cPdf_quad   = RooExponential('Pdf_quad_%s'%s_name,'Pdf_quad_%s'%s_name,rrv_x,a2)
                else:
                    N_lin       = RooRealVar('N_lin_%s'%s_name,'N_lin_%s'%s_name, (N_pos_tmp-N_neg_tmp)/2 )
                    a2_4fit     = RooRealVar('a_quad_4fit_%s'%s_name,'a_quad_4fit_%s'%s_name,-0.001,-0.01,0.)
                    a2          = RooFormulaVar('a_quad_nuis_%s'%s_name,'a_quad_nuis_%s'%s_name,'@0*@1',RooArgList(a2_4fit,self.eps))
                    a3_4fit     = RooRealVar('a_lin_4fit_%s'%s_name,'a_lin_4fit_%s'%s_name,-0.001,-0.01,0.)
                    a3          = RooFormulaVar('a_lin_nuis_%s'%s_name,'a_lin_nuis_%s'%s_name,'@0*@1',RooArgList(a3_4fit,self.eps))
                    #simple exponential sufficient for WW events
                    if   sample == 'WW':
                        cPdf_quad       = RooExponential('Pdf_quad_%s'%s_name,'Pdf_quad_%s'%s_name,rrv_x,a2)
                    elif sample == 'WZ':
                        cPdf_quad       = RooErfExpPdf('Pdf_quad_%s'%s_name,'Pdf_quad_%s'%s_name,rrv_x,a2,self.wtmp.var('Erf_offset_%s'%s_name),self.wtmp.var('Erf_width_%s'%s_name))
                a2_4fit.setConstant(kTRUE)
                a3_4fit.setConstant(kTRUE)
                #PDF for SM interference
                cPdf_lin        = RooExponential('Pdf_lin_%s'%s_name,'Pdf_lin_%s'%s_name,rrv_x,a3)

                self.Import_to_ws(self.wtmp,[cPdf_quad,cPdf_lin],1)
                self.Import_to_ws(self.wtmp,[N_quad,N_lin,scaleshape])
                
            ###make model
            #list of all coefficients
            paralist    = RooArgList(N_SM)

            #include aTGC-interference
            #read results of generator level study
            fileInter   = TFile.Open('Input/genlevel_%s_%s.root'%(sample,self.ch))
            w2          = fileInter.Get('w2')

            #get parameter values of aTGC-interference from generator level studies        
            a5_tmp      = RooRealVar('a_cwww_ccw_%s'%channel,'a_cwww_ccw_%s'%channel, w2.var('a5').getVal())
            a7_tmp      = RooRealVar('a_ccw_cb_%s'%channel,'a_ccw_cb_%s'%channel, w2.var('a7').getVal())
            a5_tmp.setConstant(kTRUE)
            a7_tmp.setConstant(kTRUE)
            #apply uncertainty parameter, bigger uncertainty for c_B in WZ
            a5          = RooFormulaVar('a_cwww_ccw_nuis_%s'%channel,'a_cwww_ccw_nuis_%s'%channel,'@0*@1',RooArgList(a5_tmp,self.eps))
            if sample=='WZ':
                a7           = RooFormulaVar('a_ccw_cb_nuis_%s'%channel,'a_ccw_cb_nuis_%s'%channel,'@0*@1',RooArgList(a7_tmp,self.eps4cbWZ))
            else:
                a7           = RooFormulaVar('a_ccw_cb_nuis_%s'%channel,'a_ccw_cb_nuis_%s'%channel,'@0*@1',RooArgList(a7_tmp,self.eps))
            
            Pdf_cwww_ccw    = RooExponential('Pdf_cwww_ccw_%s'%channel,'Pdf_cwww_ccw_%s'%channel,rrv_x,a5)
            Pdf_ccw_cb      = RooExponential('Pdf_ccw_cb_%s'%channel,'Pdf_ccw_cb_%s'%channel,rrv_x,a7)

            #get factor to scale number of events to simulation level (ratio of N_events for all atgc-parameters negative)        
            fileInHist      = TFile.Open('Output/hists4scale_%s_WV_aTGC-%s_%s.root'%(self.ch,self.binlo,self.binhi))
            hist_all3       = fileInHist.Get('c_%s_histall3'%sample)
            hist_all3.Sumw2(kTRUE)
            datahist_all3   = RooDataHist('datahist_all3','datahist_all3',RooArgList(rrv_x),hist_all3)
            getattr(self.wtmp,'import')(datahist_all3)
            fileInHist.Close()
            N_4norm         = w2.var('N_4norm').getVal()
            if options.noatgcint:
                cf = 0
            else:
                cf = datahist_all3.sumEntries() / N_4norm

            #get other coefficients
            NSM_gen     = w2.var('N_SM').getVal()
            N1220       = w2.var('N_cwww_ccw_12_20').getVal()
            N2060       = w2.var('N_ccw_cb_20_60').getVal()
            N12         = w2.var('N_cwww_12').getVal()
            N12_        = w2.var('N_cwww__12').getVal()
            N20         = w2.var('N_ccw_20').getVal()
            N20_        = w2.var('N_ccw__20').getVal()
            N60         = w2.var('N_cb_60').getVal()
            N60_        = w2.var('N_cb__60').getVal()

            ##define final coefficients, scaled by cf
            N_cwww_ccw      = RooRealVar('N_cwww_ccw_%s'%channel,'N_cwww_ccw_%s'%channel,\
                                            cf*((N1220+NSM_gen)-(N12+N20)))
            N_ccw_cb        = RooRealVar('N_ccw_cb_%s'%channel,'N_ccw_cb_%s'%channel,\
                                            cf*((N2060+NSM_gen)-(N20+N60)))

            paralist.add(RooArgList(self.wtmp.function('N_quad_%s_%s'%(self.POI[0],channel)),self.wtmp.var('cwww'),\
                                    self.wtmp.function('N_quad_%s_%s'%(self.POI[1],channel)),self.wtmp.function('N_lin_%s_%s'%(self.POI[1],channel)),self.wtmp.var('ccw'),\
                                    self.wtmp.function('N_quad_%s_%s'%(self.POI[2],channel)),self.wtmp.function('N_lin_%s_%s'%(self.POI[2],channel)),self.wtmp.var('cb')))
            paralist.add(RooArgList(N_cwww_ccw,N_ccw_cb))
            
            #parts of final signal model formula
            cwww_s      = '+@1*(@2/12)**2'
            ccw_s       = '+@3*(@5/20)**2+@4*(@5/20)'
            cb_s        = '+@6*(@8/60)**2+@7*(@8/60)'
            cwww_ccw_s  = '+@9*(@2/12)*(@5/20)'
            ccw_cb_s    = '+@10*(@5/20)*(@8/60)'
            Pdf_norm    = RooFormulaVar( 'Pdf_norm_%s'%channel, 'Pdf_norm_%s'%channel, '@0'+cwww_s+ccw_s+cb_s+cwww_ccw_s+ccw_cb_s, paralist)
            paralistN   = RooArgList()
            for i in range(11):
                paralistN.add(RooArgList(paralist.at(i)))
            paralistN.add(RooArgList(Pdf_norm))

            N1                = RooFormulaVar( 'N1_%s'%channel, 'N1_%s'%channel, '@0/@11', paralistN )
            N2                = RooFormulaVar( 'N2_%s'%channel, 'N2_%s'%channel, '(@1*(@2/12)**2)/@11', paralistN )
            #N3 ->no SM-interference for c_WWW
            N4                = RooFormulaVar( 'N4_%s'%channel, 'N4_%s'%channel, '(@3*(@5/20)**2)/@11', paralistN )
            N5                = RooFormulaVar( 'N5_%s'%channel, 'N5_%s'%channel, '(@4*(@5/20))/@11', paralistN )
            N6                = RooFormulaVar( 'N6_%s'%channel, 'N6_%s'%channel, '(@6*(@8/60)**2)/@11', paralistN )
            N7                = RooFormulaVar( 'N7_%s'%channel, 'N7_%s'%channel, '(@7*(@8/60))/@11', paralistN )
            N8                = RooFormulaVar( 'N8_%s'%channel, 'N8_%s'%channel, '(@9*(@2/12)*(@5/20))/@11', paralistN )
            #N9 ->no aTGC-interference for c_WWW/c_B
            N10                = RooFormulaVar( 'N10_%s'%channel,'N10_%s'%channel,'(@10*(@5/20)*(@8/60))/@11', paralistN )

            N_list        = RooArgList(N1,N2,N4,N5,N6,N7)
            N_list.add(RooArgList(N8,N10))
            Pdf_list        = RooArgList(SMPdf)
            Pdf_list.add(RooArgList(self.wtmp.pdf('Pdf_quad_cwww_%s'%channel),\
                                    self.wtmp.pdf('Pdf_quad_ccw_%s'%channel),self.wtmp.pdf('Pdf_lin_ccw_%s'%channel),\
                                    self.wtmp.pdf('Pdf_quad_cb_%s'%channel),self.wtmp.pdf('Pdf_lin_cb_%s'%channel)))
            Pdf_list.add(RooArgList(Pdf_cwww_ccw,Pdf_ccw_cb))
            model                = RooAddPdf('aTGC_model_%s'%channel,'aTGC_model_%s'%channel, Pdf_list, N_list)

            scale_list        = RooArgList(self.wtmp.function('scaleshape_cwww_%s'%channel), self.wtmp.function('scaleshape_ccw_%s'%channel), self.wtmp.function('scaleshape_cb_%s'%channel))
            normfactor_3d        = RooFormulaVar('normfactor_3d_%s'%channel,'normfactor_3d_%s'%channel,'1+@0+@1+@2',scale_list)
            self.wtmp.Print()

            #fit 3 pdfs for 3 atgc parameters
            for i in range(3):
                s_name        = self.POI[i] + '_' + channel
                for j in range(3):
                    self.wtmp.var(self.POI[j]).setVal(0)
                self.wtmp.var(self.POI[i]).setVal(self.PAR_MAX[self.POI[i]])

                #fit SM-interference first
                if not self.POI[i] == 'cwww' and not (sample=='WZ' and self.POI[i]=='cb'):
                    #set SM and quadratical terms to zero so only the linear term is fitted
                    N_SM_tmp = N_SM.getVal()
                    N_quad_tmp = self.wtmp.var('N_quad_%s'%s_name).getVal()
                    N_SM.setVal(0)
                    self.wtmp.var('N_quad_%s'%s_name).setVal(0)
                    
                    self.wtmp.var('a_lin_4fit_%s'%s_name).setConstant(kFALSE)
                    fitres1                = model.fitTo(self.wtmp.data('dif_datahist_%s_%s'%(sample,self.POI[i])),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer('Minuit2'))
                    self.wtmp.var('a_lin_4fit_%s'%s_name).setConstant(kTRUE)
                    self.fitresults.append(fitres1)
                    
                    N_SM.setVal(N_SM_tmp)
                    self.wtmp.var('N_quad_%s'%s_name).setVal(N_quad_tmp)

                #fit quadratical term
                self.wtmp.var('a_quad_4fit_%s'%s_name).setConstant(kFALSE)
                if sample=='WZ' and self.POI[i]!='cb':
                    self.wtmp.var('Erf_offset_%s'%s_name).setConstant(kFALSE)
                    self.wtmp.var('Erf_width_%s'%s_name).setConstant(kFALSE)
                fitres2         = model.fitTo(self.wtmp.data('pos_datahist_%s_%s'%(sample,self.POI[i])), RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
                fitres2                = model.fitTo(self.wtmp.data('pos_datahist_%s_%s'%(sample,self.POI[i])), RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer('Minuit2'))
                self.wtmp.var('a_quad_4fit_%s'%s_name).setConstant(kTRUE)
                if sample=='WZ' and self.POI[i]!='cb':
                    self.wtmp.var('Erf_offset_%s'%s_name).setConstant(kTRUE)
                    self.wtmp.var('Erf_width_%s'%s_name).setConstant(kTRUE)
                self.fitresults.append(fitres2)
                
            for i in range(3):
                self.wtmp.var(self.POI[i]).setVal(0)
            
            if options.atgc:
                #go from EFT parametrization to Lagrangian approach, taken from SI-HEP-2011-17
        
                Z_mass          = 0.0911876
                W_mass          = 0.080385
                G_F             = 11.663787
                g_weak          = math.sqrt((8*G_F*W_mass**2)/(math.sqrt(2)))
                theta_W         = math.acos(W_mass/Z_mass)
                tan_theta_W     = math.tan(theta_W)
                sin_theta_W     = math.sin(theta_W)

                coeff_cb1       = RooRealVar('coeff_cb1','coeff_cb1',2/(tan_theta_W*tan_theta_W*Z_mass*Z_mass))
                coeff_cb2       = RooRealVar('coeff_cb2','coeff_cb2',2/(sin_theta_W*sin_theta_W*Z_mass*Z_mass))
                coeff_ccw       = RooRealVar('coeff_ccw','coeff_ccw',2/(Z_mass*Z_mass))
                coeff_cwww      = RooRealVar('coeff_cwww','coeff_cwww',2/(3*g_weak*g_weak*W_mass*W_mass))

                dg1z            = RooRealVar('dg1z','dg1z',0,-1,1)
                lZ              = RooRealVar('lZ','lZ',0,-1,1)
                dkz             = RooRealVar('dkz','dkz',0,-1,1)
                dg1z.setConstant(kTRUE)
                lZ.setConstant(kTRUE)
                dkz.setConstant(kTRUE)

                cwww            = RooFormulaVar('cwww_atgc','cwww_atgc','@0*@1',RooArgList(lZ,coeff_cwww))
                ccw             = RooFormulaVar('ccw_atgc','ccw_atgc','@0*@1',RooArgList(dg1z,coeff_ccw))
                cb              = RooFormulaVar('cb_atgc','cb_atgc','@0*@1-@2*@3',RooArgList(dg1z,coeff_cb1,dkz,coeff_cb2))
                atgc_pars       = RooArgList(cwww,ccw,cb)
                N_list_atgc     = RooArgList()
                scale_list_atgc = RooArgList()

                for i in range(8):
                        customize_N        = RooCustomizer(N_list.at(i),'customize_N')
                        for j in range(3):
                                customize_N.replaceArg(self.wtmp.var(self.POI[j]),atgc_pars.at(j))
                        N_list_atgc.add(customize_N.build())
                        N_list_atgc.at(i).SetName(N_list.at(i).GetName())
                model_atgc        = RooAddPdf('aTGC_model_%s'%channel,'aTGC_model_%s'%channel, Pdf_list, N_list_atgc)
                for i in range(3):
                        customize_scale        = RooCustomizer(scale_list.at(i),'customize_scale')
                        for j in range(3):
                                customize_scale.replaceArg(self.wtmp.var(self.POI[j]),atgc_pars.at(j))
                        scale_list_atgc.add(customize_scale.build())
                        
                normfactor_3d        = RooFormulaVar('normfactor_3d_%s'%channel,'normfactor_3d_%s'%channel,'1+@0+@1+@2',scale_list_atgc)
        
                getattr(self.wtmp,'import')(model_atgc)
                getattr(self.WS,'import')(model_atgc)
                getattr(self.wtmp,'import')(normfactor_3d,RooFit.RecycleConflictNodes())
                getattr(self.WS,'import')(normfactor_3d,RooFit.RecycleConflictNodes())
                self.WS.Print()
                raw_input(channel)
            else:
                model.Print()
                self.Import_to_ws(self.wtmp,[normfactor_3d,model])
                self.Import_to_ws(self.WS,[normfactor_3d,model])
            
            #print coefficients to see contribution for all atgc-parameter positive
            for i in range(3):
                self.wtmp.var(self.POI[i]).setVal(self.PAR_MAX[self.POI[i]])
            for i in range(11):
                print paralist.at(i).GetName() + ' : ' + str(paralist.at(i).getVal())

            #print self.fitresults
            for i in range(8):
                print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())


        def Write_datacard(self):
            ### make the card for this channel and plane ID

            for region in self.regions:
                codename    = 'WWWZ_' + region + '_' + self.ch
                cfg_items   = self.cfg.items('Global')
                NlnN        = int(self.cfg.get('Global','nlnn'))
                lnN_for     = []
                lnN_name    = []
                lnN_value   = []
                bkgs_4card  = ['WJets','STop','TTbar']
                Nbkg_int    = len(bkgs_4card)
                for i in range(0,NlnN):
                    lnN_for.append([])
                    for name in self.cfg.get('Global','lnN%i_for'%(i+1)).split(','):
                        lnN_for[i].append(name)
                    lnN_name.append(self.cfg.get('Global','lnN%i_name'%(i+1)))
                    lnN_value.append([])
                    for value in self.cfg.get('Global','lnN%i_value'%(i+1)).split(','):
                        lnN_value[i].append(value)



                card = """
imax 1  number of channels
jmax {Nbkg_int}  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
            ------------""".format(Nbkg_int=Nbkg_int)
                for i in range(0,Nbkg_int):
                    card += """
shapes {bkg_name}\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS""".format(codename=codename,bkg_name=bkgs_4card[i])
                card += """
shapes data_obs\t\t\t\t\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS""".format(codename=codename)    
                card += """   
shapes ATGCPdf_{codename}\t\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS
            """.format(codename=codename)
                card += """   
------------
bin {codename} 
observation -1
------------
bin                         {codename}\t\t""".format(codename=codename)
                for i in range(0,Nbkg_int):
                    card += """\t\t\t{codename}""".format(codename=codename)
                card += """       
process\t\t\t    ATGCPdf_{codename}    """.format(codename=codename)
                for i in range(0,Nbkg_int):
                    card += """\t\t{bkg_name}""".format(bkg_name=bkgs_4card[i])
                card += """       
process                     0                        """
                for i in range(0,Nbkg_int):
                    card += """ \t\t\t\t{i}""".format(i=i+1)
                card += """       
rate                        1\t"""
                for i in range(0,Nbkg_int):
                    card += """ \t\t\t1"""
                card += """           
------------"""
                for i in range(0,NlnN):
                    # if found signal or codename in the list of names that are affected by lnN unc:
                    if (('%s_signal'%codename in lnN_for[i]) or (any(codename in s for s in lnN_for[i]))):
                        card+="""
{lnN_name} \t\t\t lnN """.format(lnN_name=lnN_name[i])
                        if ('%s_signal'%codename in lnN_for[i]):
                            # if lnN syst affects signal:
                            index_s=lnN_for[i].index('%s_signal'%codename)
                            card+=""" {lnN_value}      """.format(lnN_value=lnN_value[i][index_s])
                        else:
                            card+=""" -      """
                        for j in range(0,Nbkg_int):
                            name_for_lnN=codename
                            name_for_lnN+='_'
                            name_for_lnN+=bkgs_4card[j]
                            if (name_for_lnN in lnN_for[i]):
                                index=lnN_for[i].index(name_for_lnN)
                                card+="""\t\t\t\t{lnN_value}""".format(lnN_value=lnN_value[i][index])
                            else:
                                card+="""\t\t\t\t-"""

                print card
                cardfile = open('aC_%s.txt'%(codename),'w')
                cardfile.write(card)
                cardfile.close()
                    #cardfilenames.append('aC_%s.txt'%(codename))


        ########################
        ######MAIN CODE#########
        ########################
        def Make_input(self):

            #read data
            fileInData        = TFile.Open('Input/treeEDBR_data_xww_%s.root'%(self.ch))
            tree_data        = fileInData.Get('tree')  

            #prepare variables, parameters and temporary workspace
            if options.newtrees:
                self.Read_ATGCtree(self.ch)
            
            #make and fit signal pdf for WW and WZ
            self.Make_signal_pdf(self.rrv_mass_lvj,'WW')
            self.Make_signal_pdf(self.rrv_mass_lvj,'WZ')

            #read, rename and write bkg pdfs and bkg rates
            fileInWs        = TFile.Open('Input/wwlvj_%s_HPV_workspace.root'%self.ch)
            w_bkg        = fileInWs.Get('workspace4limit_') 

            path        ='%s/src/CombinedEWKAnalysis/CommonTools/data/anomalousCoupling'%os.environ["CMSSW_BASE"]

            for bkg in ['WJets','TTbar','STop','WW','WZ']:
                print  "rrv_number_%s_%s_mj"%(bkg,self.ch)
                getattr(self.WS,'import')(w_bkg.var("rrv_number_mj_%s_%s"%(bkg,self.ch)).clone('norm_%s_%s'%(bkg,self.ch)))

            #set shape parameter in mj-WJets floating
            w_bkg.var("rrv_c_ErfExp_WJets0_%s"%self.ch).setConstant(kFALSE)

            #import m_pruned and define ranges
            getattr(self.WS,'import')(w_bkg.var('rrv_mass_j'))
            self.WS.var('rrv_mass_j').setRange('sb_lo',40,65)
            self.WS.var('rrv_mass_j').setRange('sig',65,105)
            self.WS.var('rrv_mass_j').setRange('sb_hi',105,150)
            self.WS.var('rrv_mass_lvj').setRange(900,3500)

            #bkg-pdfs have the format '[bkg-name]_mlvj_[region]_[ch]' or '[bkg-name]_mj_[region]_[ch]'


            #create a workspace for each component in each region
            for region in self.regions:
                self.WS2 = self.WS.Clone("w")        #temporary 
                set_mj        = RooArgSet(self.WS2.var('rrv_mass_j'))
                for bkg in self.bkgs:
                    #define global norm for whole mj spectrum
                    norm_var    = RooRealVar('normvar_%s_%s'%(bkg,self.ch),'normvar_%s_%s'%(bkg,self.ch),self.WS2.var("norm_%s_%s"%(bkg,self.ch)).getVal(),0,1e4)
                    #define integral over region
                    reg_Int     = w_bkg.pdf('mj_%s_%s'%(bkg,self.ch)).createIntegral(set_mj,set_mj, region)
                    if bkg=='WJets':#norm floating for WJets, integral depends on (floating) shape parameter
                        norm_var.setConstant(kFALSE)
                        norm        = RooFormulaVar('%s_%s_%s_norm'%(bkg,region,self.ch),'%s_%s_%s_norm'%(bkg,region,self.ch),'@0*@1',RooArgList(reg_Int,norm_var))
                    else:#norm and integral fixed for rest
                        norm_var.setConstant(kTRUE)
                        norm        = RooFormulaVar('%s_%s_%s_norm'%(bkg,region,self.ch),'%s_%s_%s_norm'%(bkg,region,self.ch),'%s*@0'%reg_Int.getVal(),RooArgList(norm_var))
                    if self.regions == 'sig':
                        bkg_MWV     = w_bkg.pdf('%s_mlvj_sig_%s'%(bkg,self.ch))
                    else:
                        #pdfs from the sb fit are fitted simultaneously in the lower and upper sb
                        bkg_MWV = w_bkg.pdf('%s_mlvj_sb_%s'%(bkg,self.ch)).clone('%s_mlvj_%s_%s'%(bkg,region,self.ch))
                    #if bkg=='WJets' and 'sb' in region: #all shape parameters of the W+Jets pdf in sideband region floating ##FIXME should already be done in prepare_bkg_oneCat.py
                    #    params_mlvj   = bkg_MWV.getParameters(self.WS2.data('dataset_mj_%s'%region))
                    #    p_iter        = params_mlvj.createIterator()
                    #    p_iter.Reset()
                    #    param        = p_iter.Next()
                    #    while param:
                    #        param.setConstant(kFALSE)
                    #        param        = p_iter.Next()
                    bkg_mj          = w_bkg.pdf('%s_mj_%s_%s'%(bkg,region,self.ch))
                    #make 2d pdf
                    bkg_2d_pdf      = RooProdPdf(bkg,bkg,RooArgList(bkg_MWV,bkg_mj))
                    bkg_MWV.Print();bkg_mj.Print();bkg_2d_pdf.Print();
                    norm.SetName(bkg_2d_pdf.GetName()+'_norm')#the normalization variable must have the corresponding pdf-name + _norm
                    self.Import_to_ws(self.WS2,[bkg_2d_pdf,norm],1)

                #signal function for WW and WZ in signal region and lower/upper sideband
                #
                ##FIXME? signal function is not explicitly evaluated in the sideband region since its contribution is assumed to be negligible there
                ##FIXME need only one pdf for WW and WZ
                data_obs            = RooDataSet('data_obs','data_obs',w_bkg.data('dataset_2d_%s_%s'%(region,self.ch)),RooArgSet(self.WS2.var('rrv_mass_lvj'),self.WS2.var('mj_%s'%region)))
                getattr(self.WS2,'import')(data_obs)

                pdf_atgc_mlvj_WW    = self.WS2.pdf('aTGC_model_%s_WW'%self.ch)
                pdf_atgc_mlvj_WZ    = self.WS2.pdf('aTGC_model_%s_WZ'%self.ch)
                pdf_atgc_mj_WW      = w_bkg.pdf('m_j_WW_%s_pdf_%s'%(region,self.ch))
                pdf_atgc_mj_WZ      = w_bkg.pdf('m_j_WZ_%s_pdf_%s'%(region,self.ch))

                pdf_atgc_WW_2d  = RooProdPdf('ATGCPdf_WW_2d_%s_%s'%(region,self.ch),'ATGCPdf_WW_2d_%s_%s'%(region,self.ch),RooArgList(pdf_atgc_mlvj_WW,pdf_atgc_mj_WW))
                pdf_atgc_WZ_2d  = RooProdPdf('ATGCPdf_WZ_2d_%s_%s'%(region,self.ch),'ATGCPdf_WZ_2d_%s_%s'%(region,self.ch),RooArgList(pdf_atgc_mlvj_WZ,pdf_atgc_mj_WZ))
                #add WW and WZ pdfs to get a single signal pdf
                ##define relative coefficients for the RooAddPdf
                ###for now, sideband is scaled with scaling factor from signal region
                if True:#region == 'sig':
                    rel_norm_atgc_WW        = RooFormulaVar(pdf_atgc_WW_2d.GetName()+'_relnorm',pdf_atgc_WW_2d.GetName()+'_relnorm','@0*@1/(@0*(@1+@2))',RooArgList(self.WS2.function('normfactor_3d_%s_WW'%self.ch),self.WS2.function('WW_norm'),self.WS2.function('WZ_norm')))
                    rel_norm_atgc_WZ        = RooFormulaVar(pdf_atgc_WZ_2d.GetName()+'_relnorm',pdf_atgc_WZ_2d.GetName()+'_relnorm','@0*@2/(@0*(@1+@2))',RooArgList(self.WS2.function('normfactor_3d_%s_WZ'%self.ch),self.WS2.function('WW_norm'),self.WS2.function('WZ_norm')))
                else:
                    rel_norm_atgc_WW        = RooFormulaVar(pdf_atgc_WW_2d.GetName()+'_relnorm',pdf_atgc_WW_2d.GetName()+'_relnorm','@0/(@0+@1)',RooArgList(self.WS2.function('WW_norm'),self.WS2.function('WZ_norm')))
                    rel_norm_atgc_WZ        = RooFormulaVar(pdf_atgc_WZ_2d.GetName()+'_relnorm',pdf_atgc_WZ_2d.GetName()+'_relnorm','@1/(@0+@1)',RooArgList(self.WS2.function('WW_norm'),self.WS2.function('WZ_norm')))
                signal_2d_WWWZ      = RooAddPdf('ATGCPdf_WWWZ_%s_%s'%(region,self.ch),'ATGCPdf_WWWZ_%s_%s'%(region,self.ch),RooArgList(pdf_atgc_WW_2d,pdf_atgc_WZ_2d),RooArgList(rel_norm_atgc_WW,rel_norm_atgc_WZ))
                #final normalization
                ##FIXME aTGC not scaled in sideband region? remove it?
                signal_norm_WWWZ    = RooFormulaVar(signal_2d_WWWZ.GetName()+'_norm',signal_2d_WWWZ.GetName()+'_norm','@0*@2+@1*@3',RooArgList(self.WS2.function('normfactor_3d_%s_WW'%self.ch),self.WS2.function('normfactor_3d_%s_WZ'%self.ch),self.WS2.function('WW_norm'),self.WS2.function('WZ_norm')))

                self.Import_to_ws(self.WS2,[signal_2d_WWWZ,signal_norm_WWWZ],1)

                #output        = TFile('%s/%s_%s_%s.root'%(path,sample,region,self.ch),'recreate')
                output        = TFile('WWWZ_%s_%s_ws.root'%(region,self.ch),'recreate')
                self.WS2.SetName('proc_WWWZ_%s_%s'%(region,self.ch))
                self.WS2.Write();
                #self.WS2.pdf('aTGC_model_%s_%s'%(self.ch,sample)).Write()
                output.Close()
                print 'Write to file ' + output.GetName()


            ##create the datacards for all regions -> have to be combined with CombineCards.py
            self.Write_datacard()

            #make some plots
            if options.Make_plots:
                self.Make_plots(self.rrv_mass_lvj,'WW',self.fitresults)
                self.Make_plots(self.rrv_mass_lvj,'WZ',self.fitresults)
            
            for i in range(len(self.fitresults)):
                    self.fitresults[i].Print()

            if options.printatgc:
                self.wtmp.var('cwww').setVal(12)
                self.wtmp.var('ccw').setVal(20)
                self.wtmp.var('cb').setVal(0)
                print 'cwww and ccw positive:'
                for i in range(8):
                    print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())
                self.wtmp.var('cwww').setVal(12)
                self.wtmp.var('ccw').setVal(0)
                self.wtmp.var('cb').setVal(60)
                print 'cwww and cb positive:'
                for i in range(8):
                    print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())
                self.wtmp.var('cwww').setVal(0)
                self.wtmp.var('ccw').setVal(20)
                self.wtmp.var('cb').setVal(60)
                print 'ccw and cb positive:'
                for i in range(8):
                    print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())

                #actual yields
                for i in range(3):
                    for j in range(3):
                        self.wtmp.var(self.POI[j]).setVal(0)
                    self.wtmp.var(self.POI[i]).setVal(self.PAR_MAX[self.POI[i]])
                    print channel + ' ' + self.POI[i] + ' : ' + str(w.var('rate_VV').getVal()*normfactor_3d.getVal())

                raw_input(channel)


###run code###

if __name__ == '__main__':
    if options.chan=='elmu':
        makeWS_el        = Prepare_workspace_4limit('el',900,3500)
        makeWS_el.Make_input()
        makeWS_mu        = Prepare_workspace_4limit('mu',900,3500)
        makeWS_mu.Make_input()
    else:
        makeWS        = Prepare_workspace_4limit(options.chan,900,3500)
        gmakeWS.Make_input()
    #combine the created datacards
    output_card_name = 'aC_WWWZ_simfit'
    cmd = 'combineCards.py aC_WWWZ_sig_el.txt aC_WWWZ_sig_mu.txt aC_WWWZ_sb_lo_el.txt aC_WWWZ_sb_lo_mu.txt aC_WWWZ_sb_hi_el.txt aC_WWWZ_sb_hi_mu.txt > %s.txt'%output_card_name
    print cmd
    os.system(cmd)
    print 'generated Card : %s.txt'%output_card_name
