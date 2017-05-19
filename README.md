# aTGC2017-Signal

SIGNAL ESTIMATION  

This script creates the complete signal model and creates a datacard and a workspace per channel containing everything needed for the limit extraction.  All input files are here: `/afs/cern.ch/work/c/crenner/public/InputRootFiles/Signal`.  

use CMSSW_7_1_5

clone the repository with  
`git clone https://github.com/undeo/aTGC2017-Signal.git`  
run the script with  
`python make_PDF_input_oneCat.py -n`  
possible options:  
* `-n` read the input-trees and create RooDataHists(-> faster access); needed at the first run or when the input-trees are changed  
* `-c {channel}` only run for {channel} (el or mu)  
* `-p` make plots  
* `--savep` save the plots  
* `-b` run in batch mode  
* `--noatgcint` set atgc-interference terms to zero  
* `--printatgc` print the coefficients of the signal model  
* `--atgc` using different parametrization (lagrangian approach instead of EFT)  

the workspaces for the different channels can now be combined with  
`text2workspace.py aC_WWWZ_simfit.txt -o workspace_simfit.root -P CombinedEWKAnalysis.CommonTools.ACModel:par1par2par3_TF3_shape_Model --PO channels=WWWZ_sig_el,WWWZ_sig_mu,WWWZ_sb_lo_el,WWWZ_sb_lo_mu,WWWZ_sb_hi_el,WWWZ_sb_hi_mu --PO poi=cwww,ccw,cb --PO range_cwww=-20,20 --PO range_ccw=-30,30 --PO range_cb=-75,75`  
where  
* `-o` name of the created workspace  
* `-P` name of the used model  
* `--PO channels=` names of the channels  
* `--PO poi=` names of the paramters of interest  
* `--PO range_=` set paramter range (does't work atm but has to be added to avoid error message)  

this creates the final workspace called `workspace_simfit.root`. To inspect this workspace in ROOT you have to load `.../CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so`.

you can now run combine with (e.g. for 1d limits on cwww)  
`combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs cwww -P cwww --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-10,10 --minimizerStrategy=2 -n Example`  
where  
* `--points` is the number of scanned parameter values  
* `--redefineSignalPOIs cwww -P cwww` defines the paramter(s) of interest  
* `--freezeNuisances ccw,cb` fixes the other paramters  
* `--setPhysicsModelParameters cwww=0,ccw=0,cb=0` sets the initial parameter values  
* `--setPhysicsModelParameterRange cwww=-10,10` sets the parameter range to be scanned  
* `-n` is added to the output name  

results are saved in `higgsCombineTestExample.MultiDimFit.mH120.root`. To get 68% and 95% C.L. limits run e.g.  
`python buil1DInterval.py -10 10 higgsCombineExample.MultiDimFit.mH120.root cwww`  

to get 2dimensional limits run e.g.:  
`combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs cwww,ccw -P cwww -P ccw --freezeNuisances cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-10,10:ccw=-20,20 --minimizerStrategy=2 -n Example`  

to get the exact fit results for any point you need to run (e.g. for cwww=10)   
`combine workspace_simfit.root -M MaxLikelihoodFit --expectSignal=1 --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=10,ccw=0,cb=0 --minimizerStrategy 2 --redefineSignalPOIs cwww --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n Example`  
the output is saved in `mlfitExample.root` containing a RooFitResult `fit_s` with all final parameter values as well as a RooArgSet `norm_fit_s` with the final normalizations.  
To plot the results you can use  
`python check_combine_results.py -n Example -c mu -P cwww:10`  
which plots the mj and mlvj spectrum in all three regions while setting cwww to 10.  
If you want to get the background-only fit you have to freeze all aTGC-parameters and set a different POI, e.g.  
`combine workspace_simfit.root -M MaxLikelihoodFit --expectSignal=1 --freezeNuisances cwww,ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --minimizerStrategy 2 --redefineSignalPOIs normvar_WJets_el --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n BkgOnly`
