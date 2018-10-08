#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHybridBDTAutoPdf.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "TRandom.h"
#include "TGraph.h"
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "TFile.h"
#include "TNtuple.h"
#include "HybridGBRForest.h"
#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooDoubleCBFast.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
#include "RooLinearVar.h"
#include "RooCBExp.h"
#include "RooCBFast.h"
#include "RooGaussianFast.h"
#include "TRandom3.h"
 
using namespace RooFit;

void swapvars(RooArgList &vars) {
 
  for (int ivar=0; ivar<vars.getSize(); ++ivar) {
    TString title = vars.at(ivar)->GetTitle();
    title.ReplaceAll("ele0.","ele2.");
    title.ReplaceAll("ele1.","ele0.");
    title.ReplaceAll("ele2.","ele1.");
    vars.at(ivar)->SetTitle(title);
  }  
  
  
}


void eregtraining_data() {
  
  TString dirname = "/afs/cern.ch/work/r/rcoelhol/public/CMSSW_8_0_12/src/HiggsAnalysis/GBRLikelihood/macros/eregdata/closure/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);  
  
  //fill list of input variable expressions (any valid TTree draw expression will work here)
  std::vector<std::string> *varsf = new std::vector<std::string>;
  varsf->push_back("ele.ECALenergyRelError");
  varsf->push_back("ele.full5x5_e5x5_over_scRawEnergy");
  varsf->push_back("ele.full5x5_r9");  
  varsf->push_back("ele.full5x5_sigmaIetaIeta");
  varsf->push_back("ele.full5x5_sigmaIetaIphi");  
  varsf->push_back("ele.full5x5_sigmaIphiIphi");
  varsf->push_back("ele.full5x5_eMax_over_e5x5");
  varsf->push_back("ele.full5x5_e2nd_over_e5x5");
  varsf->push_back("ele.full5x5_eTop_over_e5x5");
  varsf->push_back("ele.full5x5_eBottom_over_e5x5");
  varsf->push_back("ele.full5x5_eLeft_over_e5x5");
  varsf->push_back("ele.full5x5_eRight_over_e5x5");
  varsf->push_back("ele.full5x5_e2x5Max_over_e5x5");
  varsf->push_back("ele.full5x5_e2x5Left_over_e5x5");
  varsf->push_back("ele.full5x5_e2x5Right_over_e5x5");
  varsf->push_back("ele.full5x5_e2x5Top_over_e5x5");
  varsf->push_back("ele.full5x5_e2x5Bottom_over_e5x5");
  varsf->push_back("ele.N_SATURATEDXTALS");
  varsf->push_back("ele.iEtaCoordinate");
  varsf->push_back("ele.iPhiCoordinate");
  varsf->push_back("ele.iEtaMod5Coordinate");
  varsf->push_back("ele.iPhiMod2Coordinate");
  varsf->push_back("ele.iEtaMod20Coordinate");
  varsf->push_back("ele.iPhiMod20Coordinate");


  std::vector<std::string> *varslist = varsf;
  
  //construct RooRealVars from variable list
  RooArgList vars_ele0;
  RooArgList vars_ele1;  
  RooArgList vars_all;  
  for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) {

    TString sname = varslist->at(ivar).c_str();
    TString title_ele0(sname);
    title_ele0.ReplaceAll("ele.","ele0.");
    TString title_ele1(sname);
    title_ele1.ReplaceAll("ele.","ele1.");

    RooRealVar *var_ele0 = new RooRealVar(TString::Format("var_ele0_%i",ivar), title_ele0 ,0.);
    RooRealVar *var_ele1 = new RooRealVar(TString::Format("var_ele1_%i",ivar), title_ele1 ,0.);
    vars_ele0.add(*var_ele0);
    vars_ele1.add(*var_ele1);
    vars_all.add(*var_ele0);
    vars_all.add(*var_ele1);      
  }
  
  vars_ele0.Print("V");
  vars_ele1.Print("V");
  vars_all.Print("V");


  //define some additional RooRealVars
  RooRealVar *tgtvar = new RooRealVar("tgtvar","mass",1.);

  const double pdgzmass = 91.1876;
  const double logpdgzmass = log(pdgzmass);
  
  RooRealVar *mass = new RooRealVar("mass","sqrt(2.0*(ele0.ECALcorrEnergy)*(ele1.ECALcorrEnergy)*(1.0-event.ECALcostheta))",90.,70.,110.);
  
  RooRealVar *logmass = new RooRealVar("logmass","0.5*log(2.0) + 0.5*log((ele0.ECALcorrEnergy)) + 0.5*log((ele1.ECALcorrEnergy)) + 0.5*log(1.0-event.ECALcostheta) -log(91.1876)",0.,log(70.)-logpdgzmass,log(110.)-logpdgzmass);
  RooRealVar *relpt_ele0 = new RooRealVar("relpt_ele0","ele0.ECALcorrPt/(sqrt(2.0*(ele0.ECALcorrEnergy)*(ele1.ECALcorrEnergy)*(1.0-event.ECALcostheta)))",1.);
  RooRealVar *relpt_ele1 = new RooRealVar("relpt_ele1","ele1.ECALcorrPt/(sqrt(2.0*(ele0.ECALcorrEnergy)*(ele1.ECALcorrEnergy)*(1.0-event.ECALcostheta)))",1.);  
  RooRealVar *eta_ele0 = new RooRealVar("eta_ele0","ele0.scEta",1.);
  RooRealVar *eta_ele1 = new RooRealVar("eta_ele1","ele1.scEta",1.);
  RooRealVar *cosdphi = new RooRealVar("cosdphi","TMath::Cos(ele0.scPhi-ele1.scPhi)",1.);
  
  RooRealVar *run = new RooRealVar("run","event.run",0.);
  RooRealVar *rhoValue = new RooRealVar("rhoValue","event.rhoValue",0.);
  
  RooRealVar *costheta = new RooRealVar("costheta","event.ECALcostheta",0.);
  RooRealVar *deltaE = new RooRealVar("deltaE","0.5*log((ele0.ECALcorrEnergy)) - 0.5*log((ele1.ECALcorrEnergy))",0.);
  RooRealVar *deltaEswap = new RooRealVar("deltaEswap","0.5*log((ele1.ECALcorrEnergy)) - 0.5*log((ele0.ECALcorrEnergy))",0.);  
  
  RooRealVar *e_ele0 = new RooRealVar("e_ele0","(ele0.ECALcorrEnergy)",0.);
  RooRealVar *scrawe_ele0 = new RooRealVar("scrawe_ele0","(ele0.scRawEnergy+ele0.scPreshowerEnergy)",0.);

  RooRealVar *e_ele1 = new RooRealVar("e_ele1","(ele1.ECALcorrEnergy)",0.);
  RooRealVar *scrawe_ele1 = new RooRealVar("scrawe_ele1","(ele1.scRawEnergy+ele1.scPreshowerEnergy)",0.);

  RooArgList vars_common;
  vars_common.add(*rhoValue);

  RooArgList vars_common_data;
  vars_common_data.add(*run);
  
  RooArgList vars_common_mass;
  vars_common_mass.add(*costheta);
  
  RooArgList varsNoEdep_ele0(vars_ele0);
  RooArgList varsNoEdep_ele1(vars_ele1);
  
  RooArgList vars;
  vars.add(vars_ele0);
  vars.add(vars_ele1);
  vars.add(vars_common);
  vars.add(vars_common_data);
  vars.add(vars_common_mass);
  vars.add(*deltaE);
  vars.add(*deltaEswap);
  vars.add(*relpt_ele0);
  vars.add(*relpt_ele1);
  vars.add(*logmass);
  vars.add(*mass);
  
  //define lists of input variables used later to define pdfs
  RooArgList corvarsNoEdep_ele0;
  corvarsNoEdep_ele0.add(varsNoEdep_ele0);
  corvarsNoEdep_ele0.add(vars_common);
  corvarsNoEdep_ele0.add(vars_common_data);
  corvarsNoEdep_ele0.add(*costheta);
  corvarsNoEdep_ele0.add(*deltaE);

  RooArgList corvarsNoEdep_ele1;
  corvarsNoEdep_ele1.add(varsNoEdep_ele1);
  corvarsNoEdep_ele1.add(vars_common);
  corvarsNoEdep_ele1.add(vars_common_data);  
  corvarsNoEdep_ele1.add(*costheta);  
  corvarsNoEdep_ele1.add(*deltaEswap);
  
  RooArgList corvars_ele0;
  corvars_ele0.add(vars_ele0);
  corvars_ele0.add(vars_common);
  corvars_ele0.add(vars_common_data);

  RooArgList corvars_ele1;
  corvars_ele1.add(vars_ele1);
  corvars_ele1.add(vars_common);
  corvars_ele1.add(vars_common_data);  
    
  RooArgList condvarsmassNoEdep;
  condvarsmassNoEdep.add(varsNoEdep_ele0);
  condvarsmassNoEdep.add(varsNoEdep_ele1);
  condvarsmassNoEdep.add(vars_common);
  condvarsmassNoEdep.add(*costheta);  
  condvarsmassNoEdep.add(*deltaE);
  
  RooArgList condvarsmassNoEdepData;
  condvarsmassNoEdepData.add(condvarsmassNoEdep);
  condvarsmassNoEdep.add(vars_common_data);    
  
      
  //define RooRealVar for event weight
  RooRealVar weightvar("weightvar","",1.);


  //load trees
  TString treeloc("een_analyzer/ElectronTree");

  TChain *tree = new TChain(treeloc);
  tree->Add("/eos/user/r/rcoelhol/80X_NTuples/DYmerged.root");

  TChain *treedata = new TChain(treeloc);
  treedata->Add("/eos/user/r/rcoelhol/80X_NTuples/datamerged.root");  

  //define cuts   
  int prescaleinit = 250;
  
  TCut selcut = "(ele0.scIsEB && ele1.scIsEB)";
  TCut selweight = "event.weight";
  
  //TCut halfweight = "(0.25*1.0)";
  
  RooDataSet *hdatanoswap = 0;
  RooDataSet *hdatanoswapD = 0;
  RooDataSet *hdataswap = 0;
  RooDataSet *hdataswapD = 0;  
  {
    //load trees into datasets using list of RooRealVars
    //RooRealVar name is arbitrary, RooRealVar Title is used as TTree draw expression (and the weight variable Title becomes the per event weight)
    weightvar.SetTitle(selcut*selweight);
    hdatanoswap = RooTreeConvert::CreateDataSet("hdatanoswap",tree,vars,weightvar);   
    hdatanoswapD = RooTreeConvert::CreateDataSet("hdatanoswapD",treedata,vars,weightvar);  
  }  
  
  //swap variable titles for electron 1 and electron 2
  swapvars(vars);

  printf("swapped vars:\n");
  for (int ivar=0; ivar<vars.getSize(); ++ivar) {
    printf("%i: %s, %s\n",ivar, vars.at(ivar)->GetName(),vars.at(ivar)->GetTitle());
  }

  
  {
    //load trees with electron 1 and 2 swapped
    printf("Loading tree");
    weightvar.SetTitle(selcut*selweight);
    hdataswap = RooTreeConvert::CreateDataSet("hdataswap",tree,vars,weightvar);   
    hdataswapD = RooTreeConvert::CreateDataSet("hdataswapD",treedata,vars,weightvar);  
  }  
   
  swapvars(vars);  

 
  
  printf("final vars:\n");
  for (int ivar=0; ivar<vars.getSize(); ++ivar) {
    printf("%i: %s, %s\n",ivar, vars.at(ivar)->GetName(),vars.at(ivar)->GetTitle());
  }
  
  
  RooArgList varsw(vars);
  varsw.add(weightvar);
  
  RooDataSet *inset = 0;
  
  
  //fill final dataset with randomized choice for "electron 0" and "electron 1"
  RooDataSet *hdata = new RooDataSet("hdata","",varsw,RooFit::WeightVar(weightvar));  
  for (int iev=0; iev<hdatanoswap->numEntries(); ++iev) {
    int dsel = RooRandom::integer(2);
    
    if (dsel==0) {
      inset = hdatanoswap;
    }
    else {
      inset = hdataswap;
    }
    
    const RooArgSet *dset = inset->get(iev);
    double weight = inset->weight();
    
    hdata->add(*dset,weight);
    
  }
  
  RooDataSet *hdataAll = hdatanoswap;
  hdataAll->append(*hdataswap);
  
  delete hdataswap;
  
  RooDataSet *hdataD = new RooDataSet("hdataD","",varsw,RooFit::WeightVar(weightvar));  
  for (int iev=0; iev<hdatanoswapD->numEntries(); ++iev) {
    int dsel = RooRandom::integer(2);
    
    if (dsel==0) {
      inset = hdatanoswapD;
    } else {
      inset = hdataswapD;
    }
    
    const RooArgSet *dset = inset->get(iev);
    double weight = inset->weight();
    
    hdataD->add(*dset,weight);
    
  }  
  
  RooDataSet *hdataAllD = hdatanoswapD;
  hdataAllD->append(*hdataswapD);  
    
  delete hdataswapD;
  
  printf("hdata:  sum = %5f, num = %i\n",hdata->sumEntries(),hdata->numEntries());
  printf("hdataD: sum = %5f, num = %i\n",hdataD->sumEntries(),hdataD->numEntries());
  
  //define  RooRealVars corresponding to non-parametric functions (scaling and smearing factors)
  RooRealVar *scalevar_ele0 = new RooRealVar("scalevar_ele0","",1.);
  RooRealVar *scalevar_ele1 = new RooRealVar("scalevar_ele1","",1.);
  RooRealVar *smearvar_ele0 = new RooRealVar("smearvar_ele0","",pow(0.005,2));
  RooRealVar *smearvar_ele1 = new RooRealVar("smearvar_ele1","",pow(0.005,2));

  RooRealVar *scaleNoEdepvar_ele0 = new RooRealVar("scaleNoEdepvar_ele0","",1.);
  RooRealVar *smearNoEdepvar_ele0 = new RooRealVar("smearNoEdepvar_ele0","",pow(0.005,2));  
  RooRealVar *scaleNoEdepvar_ele1 = new RooRealVar("scaleNoEdepvar_ele1","",1.);
  RooRealVar *smearNoEdepvar_ele1 = new RooRealVar("smearNoEdepvar_ele1","",pow(0.005,2));
  
 
  scalevar_ele0->setConstant(false);
  scalevar_ele1->setConstant(false);
  smearvar_ele0->setConstant(false);
  smearvar_ele1->setConstant(false);
  
  scaleNoEdepvar_ele0->setConstant(false);
  smearNoEdepvar_ele0->setConstant(false);
  scaleNoEdepvar_ele1->setConstant(false);
  smearNoEdepvar_ele1->setConstant(false);  
  
  //define non-parametric functions
  RooGBRFunctionFlex *scalefunc = new RooGBRFunctionFlex("scalefunc","");
  RooGBRFunctionFlex *smearfunc = new RooGBRFunctionFlex("smearfunc","");

  RooGBRFunctionFlex *scalefuncNoEdep = new RooGBRFunctionFlex("scalefuncNoEdep","");
  RooGBRFunctionFlex *smearfuncNoEdep = new RooGBRFunctionFlex("smearfuncNoEdep","");    
  
  //define mapping of input variables to non-parametric functions
  RooGBRTargetFlex *scale_ele0 = new RooGBRTargetFlex("scale_ele0","",*scalefunc,*scalevar_ele0,corvars_ele0);
  RooGBRTargetFlex *scale_ele1 = new RooGBRTargetFlex("scale_ele1","",*scalefunc,*scalevar_ele1,corvars_ele1);
  RooGBRTargetFlex *smear_ele0 = new RooGBRTargetFlex("smear_ele0","",*smearfunc,*smearvar_ele0,corvars_ele0);
  RooGBRTargetFlex *smear_ele1 = new RooGBRTargetFlex("smear_ele1","",*smearfunc,*smearvar_ele1,corvars_ele1);  

  RooGBRTargetFlex *scaleNoEdep_ele0 = new RooGBRTargetFlex("scaleNoEdep_ele0","",*scalefuncNoEdep,*scaleNoEdepvar_ele0,corvarsNoEdep_ele0);
  RooGBRTargetFlex *smearNoEdep_ele0 = new RooGBRTargetFlex("smearNoEdep_ele0","",*smearfuncNoEdep,*smearNoEdepvar_ele0,corvarsNoEdep_ele0);    
  RooGBRTargetFlex *scaleNoEdep_ele1 = new RooGBRTargetFlex("scaleNoEdep_ele1","",*scalefuncNoEdep,*scaleNoEdepvar_ele1,corvarsNoEdep_ele1);
  RooGBRTargetFlex *smearNoEdep_ele1 = new RooGBRTargetFlex("smearNoEdep_ele1","",*smearfuncNoEdep,*smearNoEdepvar_ele1,corvarsNoEdep_ele1);    
  
  const double negsmearlim = 0.0;
  
  //define transformations corresponding to parameter bounds for non-parametric outputs
  RooRealConstraint *scalelim_ele0 = new RooRealConstraint("scalelim_ele0","",*scale_ele0,0.7,1.3);
  RooRealConstraint *scalelim_ele1 = new RooRealConstraint("scalelim_ele1","",*scale_ele1,0.7,1.3);  
  RooRealConstraint *smearlim_ele0 = new RooRealConstraint("smearlim_ele0","",*smear_ele0,pow(1e-7,2),pow(0.2,2));
  RooRealConstraint *smearlim_ele1 = new RooRealConstraint("smearlim_ele1","",*smear_ele1,pow(1e-7,2),pow(0.2,2));    

  RooRealConstraint *scaleNoEdeplim_ele0 = new RooRealConstraint("scaleNoEdeplim_ele0","",*scaleNoEdep_ele0,0.7,1.3);  
  RooRealConstraint *smearNoEdeplim_ele0 = new RooRealConstraint("smearNoEdeplim_ele0","",*smearNoEdep_ele0,-pow(negsmearlim,2),pow(0.2,2));    
  RooRealConstraint *scaleNoEdeplim_ele1 = new RooRealConstraint("scaleNoEdeplim_ele1","",*scaleNoEdep_ele1,0.7,1.3);  
  RooRealConstraint *smearNoEdeplim_ele1 = new RooRealConstraint("smearNoEdeplim_ele1","",*smearNoEdep_ele1,-pow(negsmearlim,2),pow(0.2,2));     
  
  //define some additional non-parametric functions (width of pdfs for variable transformation step)
  RooRealVar *scalewidthvar_ele0 = new RooRealVar("scalewidthvar_ele0","",0.05);
  scalewidthvar_ele0->setConstant(false);
  RooGBRFunctionFlex *scalewidthfunc = new RooGBRFunctionFlex("scalewidthfunc","");
  RooGBRTargetFlex *scalewidth_ele0 = new RooGBRTargetFlex("scalewidth_ele0","",*scalewidthfunc,*scalewidthvar_ele0,corvars_ele0);
  RooRealConstraint *scalewidthlim_ele0 = new RooRealConstraint("scalewidthlim_ele0","",*scalewidth_ele0,pow(1e-7,1),pow(0.2,1));
  
  RooRealVar *smearwidthvar_ele0 = new RooRealVar("smearwidthvar_ele0","",pow(0.05,2));
  smearwidthvar_ele0->setConstant(false);
  RooGBRFunctionFlex *smearwidthfunc = new RooGBRFunctionFlex("smearwidthfunc","");
  RooGBRTargetFlex *smearwidth_ele0 = new RooGBRTargetFlex("smearwidth_ele0","",*smearwidthfunc,*smearwidthvar_ele0,corvars_ele0);
  RooRealConstraint *smearwidthlim_ele0 = new RooRealConstraint("smearwidthlim_ele0","",*smearwidth_ele0,pow(1e-4,2),pow(0.2,2));    
  
  RooRealVar *scalewidthvar_ele1 = new RooRealVar("scalewidthvar_ele1","",0.05);
  scalewidthvar_ele1->setConstant(false);
  RooGBRTargetFlex *scalewidth_ele1 = new RooGBRTargetFlex("scalewidth_ele1","",*scalewidthfunc,*scalewidthvar_ele1,corvars_ele1);
  RooRealConstraint *scalewidthlim_ele1 = new RooRealConstraint("scalewidthlim_ele1","",*scalewidth_ele1,pow(1e-7,1),pow(0.2,1));
  
  RooRealVar *smearwidthvar_ele1 = new RooRealVar("smearwidthvar_ele1","",pow(0.05,2));
  smearwidthvar_ele1->setConstant(false);
  RooGBRTargetFlex *smearwidth_ele1 = new RooGBRTargetFlex("smearwidth_ele1","",*smearwidthfunc,*smearwidthvar_ele1,corvars_ele1);
  RooRealConstraint *smearwidthlim_ele1 = new RooRealConstraint("smearwidthlim_ele1","",*smearwidth_ele1,pow(1e-4,2),pow(0.2,2));      
  
  //define non-parametric functions for ln(m) pdfs on mc and data (sums of gaussians)
  int ngaus = 6;  
  double step = 0.2/double(std::max(1,ngaus-1));
  
  RooArgList tgtsmassNoEdep;
  
  RooArgList gauspdfsNoEdep;
  RooArgList gauspdfsdataNoEdep;
  RooArgList gauspdfsdataNoEdepsingle;  
  RooArgList gauscoeffsNoEdep;
    
  for (int igaus=0; igaus<ngaus; ++igaus) {
    RooRealVar *gmeanvar = new RooRealVar(TString::Format("gNoEdep_meanvar_%i",igaus),"",log(80.)-logpdgzmass + step*igaus);
    RooRealVar *gsigmavar = new RooRealVar(TString::Format("gNoEdep_sigmavar_%i",igaus),"",0.02);
    RooRealVar *gfracvar = new RooRealVar(TString::Format("gNoEdep_fracvar_%i",igaus),"",1.0);
    
    gmeanvar->setConstant(false);
    gsigmavar->setConstant(false);
    gfracvar->setConstant(false);
    
    RooGBRFunctionFlex *gmeanfunc = new RooGBRFunctionFlex(TString::Format("gNoEdep_meanfunc_%i",igaus),"");
    RooGBRFunctionFlex *gsigmafunc = new RooGBRFunctionFlex(TString::Format("gNoEdep_sigmafunc_%i",igaus),"");
    RooGBRFunctionFlex *gfracfunc = new RooGBRFunctionFlex(TString::Format("gNoEdep_fracfunc_%i",igaus),"");
    
    RooGBRTargetFlex *gmean = new RooGBRTargetFlex(TString::Format("gNoEdep_mean_%i",igaus),"",*gmeanfunc,*gmeanvar,condvarsmassNoEdep);
    RooGBRTargetFlex *gsigma = new RooGBRTargetFlex(TString::Format("gNoEdep_sigma_%i",igaus),"",*gsigmafunc,*gsigmavar,condvarsmassNoEdep);
    RooGBRTargetFlex *gfrac = new RooGBRTargetFlex(TString::Format("gNoEdep_frac_%i",igaus),"",*gfracfunc,*gfracvar,condvarsmassNoEdep);

    RooRealConstraint *gmeanlim = new RooRealConstraint(TString::Format("gNoEdep_meanlim_%i",igaus),"",*gmean,log(75.)-logpdgzmass,log(105.)-logpdgzmass);
    RooRealConstraint *gsigmalim = new RooRealConstraint(TString::Format("gNoEdep_sigmalim_%i",igaus),"",*gsigma,sqrt(0.5*negsmearlim*negsmearlim + 1e-5),0.3);
    
    RooAbsReal *gfraclim = new RooProduct(TString::Format("gNoEdep_fraclim_%i",igaus),"",RooArgList(*gfrac,*gfrac));
     
    RooFormulaVar *gmeanscale = new RooFormulaVar(TString::Format("gNoEdep_meanscale_%i",igaus),"","@0 - 0.5*log(@1*@2)",RooArgList(*gmeanlim,*scaleNoEdeplim_ele0,*scaleNoEdeplim_ele1));
    RooFormulaVar *gsigmasmear = new RooFormulaVar(TString::Format("gNoEdep_sigmasmear_%i",igaus),"","sqrt(@0*@0 + 0.25*(@1+@2))",RooArgList(*gsigmalim,*smearNoEdeplim_ele0,*smearNoEdeplim_ele1)); 

    RooFormulaVar *gmeanscalesingle = new RooFormulaVar(TString::Format("gNoEdep_meanscalesingle_%i",igaus),"","@0 - log(@1)",RooArgList(*gmeanlim,*scaleNoEdeplim_ele0));
    RooFormulaVar *gsigmasmearsingle = new RooFormulaVar(TString::Format("gNoEdep_sigmasmearsingle_%i",igaus),"","sqrt(@0*@0 + 0.5*@1)",RooArgList(*gsigmalim,*smearNoEdeplim_ele0));     
    
    
    if (igaus==0) {
      gfraclim = new RooConstVar(TString::Format("gNoEdep_fraclimconst_%i",igaus),"",1.);
    }
    else {
      tgtsmassNoEdep.add(*gfrac);
    }
    
    RooGaussianFast *gpdf = new RooGaussianFast(TString::Format("gNoEdep_pdf_%i",igaus),"",*logmass,*gmeanlim,*gsigmalim);
    RooGaussianFast *gpdfdata = new RooGaussianFast(TString::Format("gNoEdep_pdfdata_%i",igaus),"",*logmass,*gmeanscale,*gsigmasmear);
    RooGaussianFast *gpdfdatasingle = new RooGaussianFast(TString::Format("gNoEdep_pdfdatasingle_%i",igaus),"",*logmass,*gmeanscalesingle,*gsigmasmearsingle);
    
   
    gauspdfsNoEdep.add(*gpdf);
    gauscoeffsNoEdep.add(*gfraclim);
    
    gauspdfsdataNoEdep.add(*gpdfdata);
    gauspdfsdataNoEdepsingle.add(*gpdfdatasingle);    
    
    tgtsmassNoEdep.add(*gmean);
    tgtsmassNoEdep.add(*gsigma);
    
  }  
  RooCondAddPdf masspdfNoEdep("masspdfNoEdep","",gauspdfsNoEdep,gauscoeffsNoEdep);  
  RooCondAddPdf masspdfdataNoEdep("masspdfdataNoEdep","",gauspdfsdataNoEdep,gauscoeffsNoEdep);
  RooCondAddPdf masspdfdatasingleNoEdep("masspdfdatasingleNoEdep","",gauspdfsdataNoEdepsingle,gauscoeffsNoEdep);
  

  //define lists of targets
  RooArgList tgtsNoEdep;
  tgtsNoEdep.add(tgtsmassNoEdep);
  tgtsNoEdep.add(*scaleNoEdep_ele0);
  tgtsNoEdep.add(*scaleNoEdep_ele1);
  tgtsNoEdep.add(*smearNoEdep_ele0);
  tgtsNoEdep.add(*smearNoEdep_ele1);
  
  RooArgList tgtsNoEdepsingle;
  tgtsNoEdepsingle.add(tgtsmassNoEdep);
  tgtsNoEdepsingle.add(*scaleNoEdep_ele0);
  tgtsNoEdepsingle.add(*smearNoEdep_ele0);
   
  
  RooConstVar etermconst("etermconst","",0.);  
   
  RooRealVar r("r","",1.);
  r.setConstant();


  double minweightmc = 1000.;
  double minweightdata = minweightmc*hdataD->sumEntries()/hdata->sumEntries();
  
  std::vector<double> minweights;
  minweights.push_back(100.);
  minweights.push_back(100.);
  
  //ntot.setConstant();

  TFile *fres = new TFile("fres.root","RECREATE");

  
  //run initialization fits to get starting values
  if (1) {
    std::vector<RooAbsData*> vdatainit;
    vdatainit.push_back(hdata);
    
    std::vector<RooAbsReal*> vpdfinit;
    vpdfinit.push_back(&masspdfNoEdep);
    
    RooHybridBDTAutoPdf bdtpdfinit("bdtpdfinit","",tgtsmassNoEdep,etermconst,r,vdatainit,vpdfinit);
    bdtpdfinit.SetPrescaleInit(prescaleinit);    
    bdtpdfinit.TrainForest(0);
    
  }

  if (1) {
    std::vector<RooAbsData*> vdatainit;
    vdatainit.push_back(hdata);
    vdatainit.push_back(hdataD);
    
    std::vector<RooAbsReal*> vpdfinit;
    vpdfinit.push_back(&masspdfNoEdep);
    vpdfinit.push_back(&masspdfdatasingleNoEdep);
    
    RooHybridBDTAutoPdf bdtpdfinit("bdtpdfinit2","",tgtsNoEdepsingle,etermconst,r,vdatainit,vpdfinit);
    bdtpdfinit.SetPrescaleInit(prescaleinit);    
    bdtpdfinit.TrainForest(0);
    
    scaleNoEdepvar_ele1->setVal(scaleNoEdepvar_ele0->getVal());
    scaleNoEdepvar_ele1->setError(scaleNoEdepvar_ele0->getError());
    
    smearNoEdepvar_ele1->setVal(smearNoEdepvar_ele0->getVal());
    smearNoEdepvar_ele1->setError(smearNoEdepvar_ele0->getError());      
    
  }  
  
  //run final simultaneous mass fit
  if (1) {  

    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdata);    
    vdata.push_back(hdataD);        
    
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back(&masspdfNoEdep);  
    vpdf.push_back(&masspdfdataNoEdep);         
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfNoEdep","",tgtsNoEdep,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    bdtpdfdiff.SetDoInitialFit(false);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(150);
    bdtpdfdiff.TrainForest(1500);
    
  }      
  

  //fill dataset with corrections (first step corrections which depend on dielectron correlations)
  RooFormulaVar *scaleout_ele0 = new RooFormulaVar("scaleout_ele0","","-log(@0)",RooArgList(*scaleNoEdeplim_ele0));
  
  RooRealVar *scaleoutvar_ele0 = (RooRealVar*)hdataAllD->addColumn(*scaleout_ele0);
  scaleoutvar_ele0->removeRange();
  scaleoutvar_ele0->setConstant(false);

  //define pdfs for variable transformation step
  
  RooFormulaVar *scalemean_ele0 = new RooFormulaVar("scalemean_ele0","","-log(@0)",RooArgList(*scalelim_ele0));
  RooGaussianFast *scalepdf_ele0 = new RooGaussianFast("scalepdf_ele0","",*scaleoutvar_ele0,*scalemean_ele0,*scalewidthlim_ele0);

  RooArgList tgtsscale;
  tgtsscale.add(*scale_ele0);
  tgtsscale.add(*scalewidth_ele0);
     
  
  RooRealVar *smearoutvar_ele0 = (RooRealVar*)hdataAllD->addColumn(*smearNoEdeplim_ele0);
  smearoutvar_ele0->removeRange();
  smearoutvar_ele0->setConstant(false);
  
  RooGaussianFast *smearpdf_ele0 = new RooGaussianFast("smearpdf_ele0","",*smearoutvar_ele0,*smearlim_ele0,*smearwidthlim_ele0);

  RooArgList tgtssmear;
  tgtssmear.add(*smear_ele0);
  tgtssmear.add(*smearwidth_ele0);  

  //run fit for variable transformation on scale output
  if (1) {  
    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdataAllD);    
    
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back(scalepdf_ele0);  
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff2","",tgtsscale,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    bdtpdfdiff.SetPrescaleInit(prescaleinit);    
    bdtpdfdiff.SetDoInitialFit(true);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(1500);
    bdtpdfdiff.TrainForest(1500);
  }   
  
  scale_ele1->SetUseFunc(true);
  scalewidth_ele1->SetUseFunc(true);
  

  
  
  RooWorkspace *wereg = new RooWorkspace("weregmass");
  wereg->import(*scalelim_ele0,RecycleConflictNodes());
  wereg->import(*scalelim_ele1,RecycleConflictNodes());
  wereg->import(*scalewidthlim_ele0,RecycleConflictNodes());
  wereg->import(*scalewidthlim_ele1,RecycleConflictNodes());
  wereg->import(*scaleNoEdeplim_ele0,RecycleConflictNodes());
  wereg->import(*smearNoEdeplim_ele0,RecycleConflictNodes());
  wereg->import(*scaleNoEdeplim_ele1,RecycleConflictNodes());
  wereg->import(*smearNoEdeplim_ele1,RecycleConflictNodes());  
  wereg->import(*scalepdf_ele0,RecycleConflictNodes());  
  wereg->import(masspdfNoEdep,RecycleConflictNodes());  
  wereg->import(masspdfdataNoEdep,RecycleConflictNodes());  
  
  wereg->defineSet("vars",vars,true);
  
  wereg->writeToFile("weregmass_scale.root");    
  
  //run fit for variable transformation on smearing output
  if (1) {  
    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdataAllD);    
    
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back(smearpdf_ele0);  
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff3","",tgtssmear,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    bdtpdfdiff.SetPrescaleInit(prescaleinit);    
    bdtpdfdiff.SetDoInitialFit(true);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(1500);
    bdtpdfdiff.TrainForest(1500);
  }   
  
  smear_ele1->SetUseFunc(true);
  smearwidth_ele1->SetUseFunc(true);  

  wereg->import(*smearlim_ele0,RecycleConflictNodes());  
  wereg->import(*smearlim_ele1,RecycleConflictNodes());
  wereg->import(*smearwidthlim_ele0,RecycleConflictNodes());
  wereg->import(*smearwidthlim_ele1,RecycleConflictNodes());    
  wereg->import(*smearpdf_ele0,RecycleConflictNodes());  
   
  wereg->writeToFile("weregmass.root");    

  //make some plots
  if (1) {
    TCanvas *cmc = new TCanvas;
    RooPlot *plot = logmass->frame(100);
    hdata->plotOn(plot);
    masspdfNoEdep.plotOn(plot,ProjWData(*hdata));
    plot->Draw();
    cmc->SaveAs("mplotmc.eps");

    TCanvas *cdata = new TCanvas;
    RooPlot *plotD = logmass->frame(100);
    hdataD->plotOn(plotD);
    masspdfdataNoEdep.plotOn(plotD,ProjWData(*hdataD));
    plotD->Draw();
    cdata->SaveAs("mplotdata.eps");
  }
   
  if (1) 
  {    
    TCanvas *cmod = new TCanvas;
    RooPlot *plotmod = scaleoutvar_ele0->frame(-0.02,0.02,200);
    hdataAllD->plotOn(plotmod);
    scalepdf_ele0->plotOn(plotmod,ProjWData(*hdataAllD));
    plotmod->Draw();
    cmod->SaveAs("plotmod1.eps");
  }
  
  if (1)
  {    
    TCanvas *cmod = new TCanvas;
    RooPlot *plotmod = smearoutvar_ele0->frame(0.,pow(-0.02,2),200);
    hdataAllD->plotOn(plotmod);
    smearpdf_ele0->plotOn(plotmod,ProjWData(*hdataAllD));
    plotmod->Draw();
    cmod->SaveAs("plotmod2.eps");
  }  
  
  //RooAbsArg::setDirtyInhibit(false);

  
//   RooFormulaVar *scaledmass = new RooFormulaVar("scaledmass","","sqrt(@0*@1)*@2",RooArgList(*scalelim_ele0,*scalelim_ele1,*mass));
//   RooRealVar *scaledmassvar = (RooRealVar*)hdataD->addColumn(*scaledmass);
  
//   new TCanvas;
//   RooPlot *plotm = mass->frame
  
   
  
//   RooAbsReal *condnll = masspdf.createNLL(*hdata,ConditionalObservables(condvarsmass),NumCPU(16));
//   RooAbsReal *condnllD = masspdfdata.createNLL(*hdataD,ConditionalObservables(condvars),NumCPU(16));
//   
//   printf("condnll = %5f, condnllD = %5f, totalnll = %5f\n",condnll->getVal(),condnllD->getVal(), condnll->getVal()+condnllD->getVal());
  
  
  
  //return;
  
   
   
  
  return;
  
  
}
