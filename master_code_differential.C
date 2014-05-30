#include "TStyle.h"
#include "gammaJetAnalysis/gammaJetAnalysis/commonUtility.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TProfile.h"
#include "TMath.h"

TFile* inFile1_p = 0;
TFile* inFilepp_p = 0;
TTree* inTree_p = 0;
TTree* inTreepp_p = 0;

TFile* outFile_p = 0;

const char* algType[5] = {"PuPF", "PuCalo", "VsPF", "VsCalo", "T"};                                                       

const char* fileTag1;
const char* fileTag2;

//shorthands, w/ _CFMSKIM.h                                                                                                                      
const char* finalHydjet = "finalHydjet.root";

const char* PbPbUp = "PbPb_Data_JetCutUp.root";

const char* PbPbDown = "PbPb_Data_JetCutDown.root";

const char* PbPbData = "PbPbData_new.root";

//pp stuff
const char * ppData_Down = "ppData_Down.root";

const char * Pythia80_redoPF = "Pythia80_redoPF.root";

const char * Pythia80_redo = "Pythia80_redo.root";

const char* PythiaNew = "finalPythia.root";

const char* Pythia80_QG = "pt80_pp2013_P01_prod22_v81_merged_forest_0_CFMSKIM_20140403_QG.root";

const char* PPData_UP = "PP_Data_JetCutUp.root";

const char* PPData_DOWN = "PP_Data_JetCutDown.root";

const char* Pythia80 = "Pythia80.root";

const char* DataPP= "pp_Data.root";


Float_t getDPHI( Float_t phi1, Float_t phi2) {
  Float_t dphi = phi1 - phi2;

  if ( dphi > TMath::Pi())
    dphi = dphi - 2.*(TMath::Pi());
  if ( dphi <= -(TMath::Pi()) )
    dphi = dphi + 2.*(TMath::Pi());

  if ( TMath::Abs(dphi) > TMath::Pi()) {
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than TMath::Pi() " << std::endl;
  }

  return dphi;
}

Float_t getAbsDphi( Float_t phi1, Float_t phi2) {
  return TMath::Abs( getDPHI(phi1, phi2) ) ;
}

Bool_t sameSign(Double_t num1, Double_t num2){
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}

void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY)
{
  handsomeTH1N(uglyTH1);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->SetNdivisions(ndivX);
  uglyTH1->SetNdivisions(ndivY, "Y");
}

void niceTProf(TProfile* uglyTProf, float max , float min, float ndivX, float ndivY)
{
  uglyTProf->GetYaxis()->SetTitleOffset(1.25);
  uglyTProf->GetXaxis()->CenterTitle();
  uglyTProf->GetYaxis()->CenterTitle();
  uglyTProf->SetMaximum(max);
  uglyTProf->SetMinimum(min);
  uglyTProf->SetNdivisions(ndivX);
  uglyTProf->SetNdivisions(ndivY, "Y");

  uglyTProf->SetMarkerColor(1);
  uglyTProf->SetMarkerSize(1);
  uglyTProf->SetMarkerStyle(20);
  uglyTProf->SetLineColor(1);
}

void niceTGraphErrors(TGraphErrors* uglyTGraph, float max, float min)
{
  uglyTGraph->GetYaxis()->SetTitleOffset(1.25);
  uglyTGraph->GetXaxis()->CenterTitle();
  uglyTGraph->GetYaxis()->CenterTitle();
  uglyTGraph->SetMaximum(max);
  uglyTGraph->SetMinimum(min);

  uglyTGraph->SetMarkerColor(1);
  uglyTGraph->SetMarkerSize(1);
  uglyTGraph->SetMarkerStyle(20);
  uglyTGraph->SetLineColor(1);
}

// 0 == PuPF, 1 == PuCalo, 2 == VsPF, 3 == VsCalo, 4 == Truth

TCut makeSetCut(Int_t setNum)
{
  TCut setCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeSetCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return setCut;
  }

  setCut = Form("eventSet[%d]", setNum);

  return setCut;
}


TCut makeCentCut(Int_t centLow, Int_t centHi)
{
  TCut centCut = "";
  if(centLow >= 0 && centHi >= centLow && centHi <= 199)
    centCut = Form("hiBin >= %d && hiBin <= %d", centLow, centHi);
  else
    std::cout << "makeCentCut: centLow/centHi incorrectly specified, empty cut returned" << std::endl;

  return centCut;
}


TCut makeAsymmCut(Int_t setNum, Float_t asymmLow, Float_t asymmHi)
{
  TCut asymmCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeAsymmCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return asymmCut;
  }

  if(asymmLow >= .00 && asymmHi >= asymmLow && asymmHi <= 1.)
    asymmCut = Form("AlgJtAsymm[%d] > %f && AlgJtAsymm[%d] < %f ", setNum, asymmLow, setNum, asymmHi);
  else
    std::cout << "makeAsymmCut: asymmLow/asymmHi incorrectly specified, empty cut returned" << std::endl;

  return asymmCut;
}


TCut makeEtaCut(Int_t setNum, Float_t overallCut = 2.0, const char* GLN = "N")
{
  TCut etaCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeEtaCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return etaCut;
  }

  const char* leadJt = Form("AlgLeadJtEta[%d]", setNum);
  const char* subLeadJt = Form("AlgSubLeadJtEta[%d]", setNum);

  etaCut = Form("TMath::Abs(%s) < %f && TMath::Abs(%s) < %f", leadJt, overallCut, subLeadJt, overallCut);

  if(strcmp(GLN, "N") == 0)
    return etaCut;

  if(strcmp(GLN, "G") == 0){
    etaCut = Form("(TMath::Abs(%s) > 1.0 || TMath::Abs(%s) > 1.0) && TMath::Abs(%s) < %f && TMath::Abs(%s) < %f", leadJt, subLeadJt, leadJt, overallCut, subLeadJt, overallCut);
  }
  else if(strcmp(GLN, "L") == 0){
    etaCut = Form("TMath::Abs(%s) < 1.0 && TMath::Abs(%s) < 1.0", leadJt, subLeadJt);
  }
  else
    std::cout << "makeEtaCut: GLN incorrectly specified, returning blank cut" << std::endl;

  return etaCut;
}



TCut makeDelPhiCut(Int_t setNum, Float_t delPhiLow = 2*TMath::Pi()/3)
{
  TCut delPhiCut = "";

  if(setNum > 4 || setNum < 0){
    std::cout << "makeDelPhiCut: setNum must be between 0-4, empty cut returned" << std::endl;
    return delPhiCut;
  }

  const char* jtDelPhi = Form("AlgJtDelPhi[%d]", setNum);

  delPhiCut = Form("%s > %f", jtDelPhi, delPhiLow);

  return delPhiCut;
}

void makeImbDelRGraphPP(TTree* getTree_p, const char* outName, const char* gorr, Int_t setNum, const char* perpProj, const char* FPT, Int_t Aj_low, Int_t Aj_high, Int_t graphLow, Int_t graphHi, const char* GLN = "N", const char* Corr = "")
{
  inFile1_p->cd();

  Int_t setCorrNum = setNum;
  if(strcmp("", Corr) != 0)
    setCorrNum = setNum + 5;

  const char* title = Form("%s%sImbDelR%s%s%s_%s_Aj%d%d_PP", gorr, algType[setNum], perpProj, FPT, Corr, GLN, (Int_t)(Aj_low), (Int_t)(Aj_high));

  TGraphErrors* imbDelRGraph_p = new TGraphErrors(10);
  imbDelRGraph_p->GetXaxis()->SetLimits(0.00, 2.00);
  niceTGraphErrors(imbDelRGraph_p, graphHi, graphLow);

  TH1F* getHist_p;

  const char* coneR[10] = {"1C", "2C", "3C", "4C", "5C", "6C", "7C", "8C", "9C", "10C"};

  TCut setCut = makeSetCut(setNum);
  //TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 0.5, GLN);
  TCut phiCut = makeDelPhiCut(setNum,5*TMath::Pi()/6.0);
  TCut assymCut = makeAsymmCut(setNum,Aj_low/100.0,Aj_high/100.0);
  TCut jetLCut = Form("AlgLeadJtPt[%d] > 120", setNum);
  TCut jetsubLCut = Form("AlgSubLeadJtPt[%d] > 50", setNum);
  //TCut gluonCut = Form("isGluonJet[%d]",setNum);
  //TCut quarkCut = Form("isQuarkJet[%d]",setNum);

  const char* name1[10] = {"1C(10000, -10000, 10000)", "2C(10000, -10000, 10000)", "3C(10000, -10000, 10000)", "4C(10000, -10000, 10000)", "5C(10000, -10000, 10000)", "6C(10000, -10000, 10000)", "7C(10000, -10000, 10000)", "8C(10000, -10000, 10000)", "9C(10000, -10000, 10000)", "10C(10000, -10000, 10000)"};

  Float_t point[10] = {.1, .3, .5, .7, .9, 1.1, 1.3, 1.5, 1.7, 1.9};
  Float_t xErr[10] = {.1, .1, .1, .1, .1, .1, .1, .1, .1, .1};
  Float_t netimb=0.0;
  Float_t netimberr=0.0;

  for(Int_t binIter = 0; binIter < 10; binIter++){
    TString var = Form("%sAlgImb%s%s%s[%d]", gorr, perpProj, coneR[binIter], FPT, setCorrNum);
    getTree_p->Project(name1[binIter], var, setCut && etaCut && phiCut && assymCut && jetLCut && jetsubLCut);

    getHist_p = (TH1F*)inFile1_p->Get(coneR[binIter]);
      imbDelRGraph_p->SetPoint(binIter, point[binIter], getHist_p->GetMean());
      imbDelRGraph_p->SetPointError(binIter, xErr[binIter], getHist_p->GetMeanError());
  }

  imbDelRGraph_p->GetXaxis()->SetLimits(0.00, 2.00);
  niceTGraphErrors(imbDelRGraph_p, graphHi, graphLow);

  outFile_p = new TFile(outName, "UPDATE");
  imbDelRGraph_p->Write(title);
  outFile_p->Close();

  delete outFile_p;
  delete imbDelRGraph_p;
}

void makeImbDelRGraph(TTree* getTree_p, const char* outName, const char* gorr, Int_t setNum, const char* perpProj, const char* FPT, Int_t centLow, Int_t centHi, Int_t Aj_low, Int_t Aj_high, Int_t graphLow, Int_t graphHi, const char* GLN = "N", const char* Corr = "")
{
  inFile1_p->cd();

  Int_t setCorrNum = setNum;
  if(strcmp("", Corr) != 0)
    setCorrNum = setNum + 5;

  const char* title = Form("%s%sImbDelR%s%s%s_%d%d_%s_%s_Aj%d%d_g", gorr, algType[setNum], perpProj, FPT, Corr, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), GLN, fileTag1, (Int_t)(Aj_low), (Int_t)(Aj_high));

  TGraphErrors* imbDelRGraph_p = new TGraphErrors(10);
  imbDelRGraph_p->GetXaxis()->SetLimits(0.00, 2.00);
  niceTGraphErrors(imbDelRGraph_p, graphHi, graphLow);

  TH1F* getHist_p;

  const char* coneR[10] = {"1C", "2C", "3C", "4C", "5C", "6C", "7C", "8C", "9C", "10C"};

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 0.5, GLN);
  TCut phiCut = makeDelPhiCut(setNum,5*TMath::Pi()/6.0);
  TCut assymCut = makeAsymmCut(setNum,Aj_low/100.0,Aj_high/100.0);
  //TCut gluonCut = Form("isGluonJet[%d]",setNum);
  //TCut quarkCut = Form("isQuarkJet[%d]",setNum);
  TCut jetLCut = Form("AlgLeadJtPt[%d] > 120", setNum);
  //TCut pthatCut = "pthat>80";

  const char* name1[10] = {"1C(10000, -10000, 10000)", "2C(10000, -10000, 10000)", "3C(10000, -10000, 10000)", "4C(10000, -10000, 10000)", "5C(10000, -10000, 10000)", "6C(10000, -10000, 10000)", "7C(10000, -10000, 10000)", "8C(10000, -10000, 10000)", "9C(10000, -10000, 10000)", "10C(10000, -10000, 10000)"};

  Float_t point[10] = {.1, .3, .5, .7, .9, 1.1, 1.3, 1.5, 1.7, 1.9};
  Float_t xErr[10] = {.1, .1, .1, .1, .1, .1, .1, .1, .1, .1};
  Float_t netimb=0.0; 
  Float_t netimberr=0.0; 

  for(Int_t binIter = 0; binIter < 10; binIter++){
    TString var = Form("%sAlgImb%s%s%s[%d]", gorr, perpProj, coneR[binIter], FPT, setCorrNum);
    //getTree_p->Project(name1[binIter], var, Form("pthatWeight*centWeight_hatAll_0[3]")*(setCut && centCut && etaCut && phiCut && jetLCut && assymCut && pthatCut && gluonCut));
      getTree_p->Project(name1[binIter], var,(setCut && centCut && etaCut && phiCut && jetLCut && assymCut));

    getHist_p = (TH1F*)inFile1_p->Get(coneR[binIter]);
      imbDelRGraph_p->SetPoint(binIter, point[binIter], getHist_p->GetMean());
      imbDelRGraph_p->SetPointError(binIter, xErr[binIter], getHist_p->GetMeanError());
  }

  imbDelRGraph_p->GetXaxis()->SetLimits(0.00, 2.00);
  niceTGraphErrors(imbDelRGraph_p, graphHi, graphLow);

  outFile_p = new TFile(outName, "UPDATE");
  imbDelRGraph_p->Write(title);

  outFile_p->Close();

  delete outFile_p;
  delete imbDelRGraph_p;
}

Double_t sumYForPTStack(Double_t dIn = 0, Double_t comp1 = 0, Double_t comp2 = 0, Double_t comp3 = 0, Double_t comp4 = 0)
{
  Double_t dOut = dIn;

  if(sameSign(comp1, dOut))
    dOut += comp1;

  if(sameSign(comp2, dOut))
    dOut += comp2;

  if(sameSign(comp3, dOut))
    dOut += comp3;

  if(sameSign(comp4, dOut))
    dOut += comp4;

  return dOut;
}

void makeHistForPtStack(TGraph* g0_1_p, TGraph* g1_2_p, TGraph* g2_4_p, TGraph* g4_8_p, TGraph* g8_100_p, TGraph* gF_p, TH1F* h0_1_p, TH1F* h1_2_p, TH1F* h2_4_p, TH1F* h4_8_p, TH1F* h8_100_p, TH1F* hF_p, Int_t pos = 4)
{
  Int_t points = 10;

  Double_t x0_1[points];
  Double_t y0_1[points];
  Double_t x1_2[points];
  Double_t y1_2[points];
  Double_t x2_4[points];
  Double_t y2_4[points];
  Double_t x4_8[points];
  Double_t y4_8[points];
  Double_t x8_100[points];
  Double_t y8_100[points];
  Double_t xF[points];
  Double_t yF[points];

  for(Int_t iter = 0; iter < points; iter++){

    g0_1_p->GetPoint(iter, x0_1[iter], y0_1[iter]);
    g1_2_p->GetPoint(iter, x1_2[iter], y1_2[iter]);
    g2_4_p->GetPoint(iter, x2_4[iter], y2_4[iter]);
    g4_8_p->GetPoint(iter, x4_8[iter], y4_8[iter]);
    g8_100_p->GetPoint(iter, x8_100[iter], y8_100[iter]);
    gF_p->GetPoint(iter, xF[iter], yF[iter]);

    h8_100_p->SetBinContent(iter + 1, y8_100[iter]);
    h8_100_p->SetBinError(iter + 1, g8_100_p->GetErrorY(iter));

    h4_8_p->SetBinContent(iter + 1, sumYForPTStack(y4_8[iter], y8_100[iter]));
    h4_8_p->SetBinError(iter + 1, g4_8_p->GetErrorY(iter));

    h2_4_p->SetBinContent(iter + 1, sumYForPTStack(y2_4[iter], y4_8[iter], y8_100[iter]));
    h2_4_p->SetBinError(iter + 1, g2_4_p->GetErrorY(iter));

    h1_2_p->SetBinContent(iter + 1, sumYForPTStack(y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h1_2_p->SetBinError(iter + 1, g1_2_p->GetErrorY(iter));

    h0_1_p->SetBinContent(iter + 1, sumYForPTStack(y0_1[iter]+y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h0_1_p->SetBinError(iter + 1, g0_1_p->GetErrorY(iter));

    hF_p->SetBinContent(iter + 1, yF[iter]);
    hF_p->SetBinError(iter + 1, gF_p->GetErrorY(iter));

    if(iter == 7){
      std::cout << std::endl;

      std::cout << h8_100_p->GetBinContent(iter+1) << std::endl;
      std::cout << h4_8_p->GetBinContent(iter+1) << std::endl;
      std::cout << h2_4_p->GetBinContent(iter+1) << std::endl;
      std::cout << h1_2_p->GetBinContent(iter+1) << std::endl;
      std::cout << h0_1_p->GetBinContent(iter+1) << std::endl;
      std::cout << hF_p->GetBinContent(iter+1) << std::endl;

      std::cout << std::endl;
    }
  }

  h8_100_p->SetXTitle(Form("#Delta R"));
  h4_8_p->SetXTitle(Form("#Delta R"));
  h2_4_p->SetXTitle(Form("#Delta R"));
  h1_2_p->SetXTitle(Form("#Delta R"));
  h0_1_p->SetXTitle(Form("#Delta R"));
  hF_p->SetXTitle(Form("#Delta R"));

  if(pos == 1){
    h8_100_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h4_8_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h2_4_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h1_2_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    h0_1_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
    hF_p->SetYTitle("<#slash{p}_{T}^{||}> (GeV/c)");
  }
  return;
}

Double_t quadSum(Double_t one, Double_t two)
{
  Double_t err = TMath::Sqrt(one*one + two*two);
  return err;
}

void makeSubHistForPtStack(TGraph* g0_1_p, TGraph* g1_2_p, TGraph* g2_4_p, TGraph* g4_8_p, TGraph* g8_100_p, TGraph* gF_p, TH1F* h0_1_p, TH1F* h1_2_p, TH1F* h2_4_p, TH1F* h4_8_p, TH1F* h8_100_p, TH1F* hF_p, TGraph* pp0_1_p, TGraph* pp1_2_p, TGraph* pp2_4_p, TGraph* pp4_8_p, TGraph* pp8_100_p, TGraph* ppF_p, Int_t pos = 4, Bool_t disp = false)
{
  Int_t points = 10;

  Double_t x0_1[points];
  Double_t y0_1[points];
  Double_t x1_2[points];
  Double_t y1_2[points];
  Double_t x2_4[points];
  Double_t y2_4[points];
  Double_t x4_8[points];
  Double_t y4_8[points];
  Double_t x8_100[points];
  Double_t y8_100[points];
  Double_t xF[points];
  Double_t yF[points];

  Double_t xPP0_1[points];
  Double_t yPP0_1[points];
  Double_t xPP1_2[points];
  Double_t yPP1_2[points];
  Double_t xPP2_4[points];
  Double_t yPP2_4[points];
  Double_t xPP4_8[points];
  Double_t yPP4_8[points];
  Double_t xPP8_100[points];
  Double_t yPP8_100[points];
  Double_t xPPF[points];
  Double_t yPPF[points];

  for(Int_t iter = 0; iter < points; iter++){
    g0_1_p->GetPoint(iter, x0_1[iter], y0_1[iter]);
    g1_2_p->GetPoint(iter, x1_2[iter], y1_2[iter]);
    g2_4_p->GetPoint(iter, x2_4[iter], y2_4[iter]);
    g4_8_p->GetPoint(iter, x4_8[iter], y4_8[iter]);
    g8_100_p->GetPoint(iter, x8_100[iter], y8_100[iter]);
    gF_p->GetPoint(iter, xF[iter], yF[iter]);

    pp0_1_p->GetPoint(iter, xPP0_1[iter], yPP0_1[iter]);
    pp1_2_p->GetPoint(iter, xPP1_2[iter], yPP1_2[iter]);
    pp2_4_p->GetPoint(iter, xPP2_4[iter], yPP2_4[iter]);
    pp4_8_p->GetPoint(iter, xPP4_8[iter], yPP4_8[iter]);
    pp8_100_p->GetPoint(iter, xPP8_100[iter], yPP8_100[iter]);
    ppF_p->GetPoint(iter, xPPF[iter], yPPF[iter]);

    if(iter == 3 && disp){
      std::cout << yF[iter] << ", " << yPPF[iter] << ", " << yF[iter] - yPPF[iter] << std::endl;
    }


    y0_1[iter] += -yPP0_1[iter];
    y1_2[iter] += -yPP1_2[iter];
    y2_4[iter] += -yPP2_4[iter];
    y4_8[iter] += -yPP4_8[iter];
    y8_100[iter] += -yPP8_100[iter];
    yF[iter] += -yPPF[iter];

    if(iter == 3)
      std::cout << yF[iter] << std::endl;


    h8_100_p->SetBinContent(iter + 1, y8_100[iter]);
    h8_100_p->SetBinError(iter + 1, quadSum(g8_100_p->GetErrorY(iter), pp8_100_p->GetErrorY(iter)));
    
    h4_8_p->SetBinContent(iter + 1, sumYForPTStack(y4_8[iter], y8_100[iter]));
    h4_8_p->SetBinError(iter + 1, quadSum(g4_8_p->GetErrorY(iter), pp4_8_p->GetErrorY(iter)));
    
    h2_4_p->SetBinContent(iter + 1, sumYForPTStack(y2_4[iter], y4_8[iter], y8_100[iter]));
    h2_4_p->SetBinError(iter + 1, quadSum(g2_4_p->GetErrorY(iter), pp2_4_p->GetErrorY(iter)));

    h1_2_p->SetBinContent(iter + 1, sumYForPTStack(y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h1_2_p->SetBinError(iter + 1, quadSum(g1_2_p->GetErrorY(iter), pp1_2_p->GetErrorY(iter)));

    h0_1_p->SetBinContent(iter + 1, sumYForPTStack(y0_1[iter], y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
    h0_1_p->SetBinError(iter + 1, quadSum(g0_1_p->GetErrorY(iter), pp0_1_p->GetErrorY(iter)));

    hF_p->SetBinContent(iter + 1, yF[iter]);
    hF_p->SetBinError(iter + 1, quadSum(gF_p->GetErrorY(iter), ppF_p->GetErrorY(iter)));

  }
  return;
}

void drawHistToPTStack(TH1F* drawHist_p, Int_t color, const char* drawOpt)
{
  drawHist_p->SetFillColor(color);
  drawHist_p->SetMarkerStyle(6);
  drawHist_p->SetMarkerSize(.5);
  drawHist_p->Draw(drawOpt);
  drawHist_p->DrawCopy("E1 SAME");
}

void makeImbDelRPtStack(const char* fileName, const char* gorr, Int_t setNum, Int_t Aj_low, Int_t Aj_high, const char* perpProj, const char* GLN, const char* Corr = "")
{
  TFile* panelFile_p = new TFile(fileName, "UPDATE");

  TGraphErrors* getGraph1_p[6];
  TGraphErrors* getGraph2_p[6];
  TGraphErrors* getGraph_pp_p[6];

  Float_t binArrayX[11] = {.00, .20, .40, .60, .80, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

  TH1F* hist1_p[6];
  TH1F* hist2_p[6];
  TH1F* hist1_sub_p[6];
  TH1F* hist2_sub_p[6];
  TH1F* hist_pp_p[6];

  TH1F* dummy_p = new TH1F("dummy", "dummy", 10, binArrayX);
  niceTH1(dummy_p, 7.5, -7.5, 505, 406);

  for(Int_t iter = 0; iter < 10; iter++){
    dummy_p->SetBinContent(iter + 1, 10);
  }

  dummy_p->SetYTitle("(PYTHIA+HYDJET)-PYTHIA");

  const char* namePT[6] = {"0_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 0; histIter < 6; histIter++){
    hist1_p[histIter] = new TH1F(Form("hist1%s_Aj%d%d_p", namePT[histIter], Aj_low, Aj_high), Form("hist1%s_Aj%d%d_p", namePT[histIter], Aj_low, Aj_high), 10, binArrayX);
    hist2_p[histIter] = new TH1F(Form("hist2%s_Aj%d%d_p", namePT[histIter], Aj_low, Aj_high), Form("hist2%s_Aj%d%d_p", namePT[histIter], Aj_low, Aj_high), 10, binArrayX);
    hist_pp_p[histIter] = new TH1F(Form("hist_pp%s_Aj%d%d_p", namePT[histIter], Aj_low, Aj_high), Form("hist_pp%s_Aj%d%d_p", namePT[histIter], Aj_low, Aj_high), 10, binArrayX);
    hist1_sub_p[histIter] = new TH1F(Form("hist1_sub_%s_%s_%d_Aj%d%dp",gorr, namePT[histIter],setNum, Aj_low, Aj_high), Form("hist1_sub_%s_%s_%d_Aj%d%dp", gorr, namePT[histIter],setNum, Aj_low, Aj_high), 10, binArrayX);
    hist2_sub_p[histIter] = new TH1F(Form("hist2_sub_%s_%s_%d_Aj%d%dp",gorr, namePT[histIter],setNum, Aj_low, Aj_high), Form("hist2_sub_%s_%s_%d_Aj%d%dp", gorr, namePT[histIter],setNum, Aj_low, Aj_high), 10, binArrayX);

    niceTH1(hist1_p[histIter], 50, -50, 505, 406);
    niceTH1(hist2_p[histIter], 50, -50, 505, 406);
    niceTH1(hist_pp_p[histIter], 50, -50, 505, 406);
    niceTH1(hist2_sub_p[histIter], 7.5, -7.5, 505, 406);
    niceTH1(hist1_sub_p[histIter], 7.5, -7.5, 505, 406);

    getGraph1_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbDelR%s%s%s_30100_%s_%s_Aj%d%d_g", gorr, algType[setNum], perpProj, namePT[histIter], Corr, GLN, fileTag1, Aj_low, Aj_high));
    getGraph2_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbDelR%s%s%s_030_%s_%s_Aj%d%d_g", gorr, algType[setNum], perpProj, namePT[histIter], Corr, GLN, fileTag1, Aj_low, Aj_high));
    getGraph_pp_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbDelR%s%s%s_%s_Aj%d%d_PP", gorr, algType[setNum], perpProj, namePT[histIter], Corr, GLN, Aj_low, Aj_high));
  }

  makeHistForPtStack(getGraph_pp_p[0], getGraph_pp_p[1], getGraph_pp_p[2], getGraph_pp_p[3], getGraph_pp_p[4], getGraph_pp_p[5], hist_pp_p[0], hist_pp_p[1], hist_pp_p[2], hist_pp_p[3], hist_pp_p[4], hist_pp_p[5], 1);
  makeSubHistForPtStack(getGraph1_p[0], getGraph1_p[1], getGraph1_p[2], getGraph1_p[3], getGraph1_p[4], getGraph1_p[5], hist1_sub_p[0], hist1_sub_p[1], hist1_sub_p[2], hist1_sub_p[3], hist1_sub_p[4], hist1_sub_p[5], getGraph_pp_p[0], getGraph_pp_p[1], getGraph_pp_p[2], getGraph_pp_p[3], getGraph_pp_p[4], getGraph_pp_p[5]);
  makeSubHistForPtStack(getGraph2_p[0], getGraph2_p[1], getGraph2_p[2], getGraph2_p[3], getGraph2_p[4], getGraph2_p[5], hist2_sub_p[0], hist2_sub_p[1], hist2_sub_p[2], hist2_sub_p[3], hist2_sub_p[4], hist2_sub_p[5], getGraph_pp_p[0], getGraph_pp_p[1], getGraph_pp_p[2], getGraph_pp_p[3], getGraph_pp_p[4], getGraph_pp_p[5]);   


  TCanvas* profPanel_p;
  profPanel_p = new TCanvas(Form("%s%sImbDelR%s%sPTStack_%s_%s_Aj%d%d_c", gorr, algType[setNum], perpProj, Corr, GLN, fileTag1, Aj_low, Aj_high), Form("%s%sImbDelR%s%sPTStack_%s_%s_Aj%d%d_c", gorr, algType[setNum], perpProj, Corr, GLN, fileTag1, Aj_low, Aj_high), 800, 600);
  profPanel_p->Divide(3, 2, 0, 0);

  TLegend* leg;
  if(strcmp(gorr, "g") == 0)
   // leg = new TLegend(0.3, 0.55, 0.97, 0.95, "Truth #slash{p}_{T}^{||}, A_{J}}<0.22");
  //else
   // leg = new TLegend(0.3, 0.55, 0.97, 0.95, Form("%sTrk #slash{p}_{T}^{||}, A_{J}<0.22", Corr));


leg = new TLegend(0.3, 0.55, 0.97, 0.95, "Truth #slash{p}_{T}^{||}");
  else
    leg = new TLegend(0.3, 0.55, 0.97, 0.95, Form("%sTrk #slash{p}_{T}^{||}", Corr));
  //leg = new TLegend(0.15, 0.6, 0.95, 0.99, Form("Truth #slash{p}_{T}^{||} v. %s #Delta R, %s", algType[setNum], fileTag1));
  //else
    //leg = new TLegend(0.15, 0.6, 0.95, 0.95, Form("%sTrk #slash{p}_{T}^{||} v. %s #Delta R, %s", Corr, algType[setNum], fileTag1));


  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(.05);
  leg->SetBorderSize(0);
  leg->SetEntrySeparation(0.08);

  profPanel_p->cd(1);
  drawHistToPTStack(hist_pp_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(hist_pp_p[1], kYellow - 9, "E1 HIST SAME");  
  drawHistToPTStack(hist_pp_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(hist_pp_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(hist_pp_p[4], kRed + 1, "E1 HIST SAME");

  hist_pp_p[5]->SetMarkerStyle(25);
  hist_pp_p[5]->SetMarkerColor(4);
  hist_pp_p[5]->DrawCopy("SAME E1");

  profPanel_p->cd(4);

  leg->AddEntry(hist_pp_p[0], ".5 < p_{T} < 1", "F");
  leg->AddEntry(hist_pp_p[1], "1 < p_{T} < 2", "F");
  leg->AddEntry(hist_pp_p[2], "2 < p_{T} < 4", "F");
  leg->AddEntry(hist_pp_p[3], "4 < p_{T} < 8", "F");
  leg->AddEntry(hist_pp_p[4], "8 < p_{T}", "F");
  leg->AddEntry(hist_pp_p[5], "integrated #slash{p}_{T,pp}^{||}", "P");

  TLegend *  leg2 = new TLegend(.15,0.25,.9,.45,"");
  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.05);
  leg2->SetBorderSize(0);
  leg->SetEntrySeparation(0.05);
  leg2->AddEntry((TObject*)0,"Lead Jet p_{T} > 120 GeV/c","");
  leg2->AddEntry((TObject*)0,"Sublead Jet p_{T} > 50 GeV/c","");
  leg2->AddEntry((TObject*)0,"Lead/Sublead Jet |#eta| < .5","");
  leg2->AddEntry((TObject*)0,"Jet #Delta #phi > 5#pi/6","");
  leg2->AddEntry((TObject*)0,"Gluon leading jets","");

  dummy_p->Draw("HIST");

  leg->Draw("SAME");
  leg2->Draw("SAME");

  profPanel_p->cd(1);

  TLine* zeroLine_p = new TLine(0., 0., 2.0, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(1);
  zeroLine_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetTextSize(0.06);
  label_p->SetNDC();
  label_p->DrawLatex(.2,.8,"CMS Preliminary");
  //label_p->DrawLatex(.2,.87,"pp, #sqrt{s_{NN}} = 2.76 TeV, 5.3/pb");
  label_p->DrawLatex(0.2,0.87,"PYTHIA");
  
  profPanel_p->cd(2);

  makeHistForPtStack(getGraph1_p[0], getGraph1_p[1], getGraph1_p[2], getGraph1_p[3], getGraph1_p[4], getGraph1_p[5], hist1_p[0], hist1_p[1], hist1_p[2], hist1_p[3], hist1_p[4], hist1_p[5], 2);

  drawHistToPTStack(hist1_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(hist1_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(hist1_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(hist1_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(hist1_p[4], kRed + 1, "E1 HIST SAME");

  hist1_p[5]->DrawCopy("SAME E1");

  zeroLine_p->Draw();

  label_p->DrawLatex(.2, .80, "30-100%");
  //label_p->DrawLatex(.1,.87,"PbPb, #sqrt{s_{NN}} = 2.76 TeV, 150/#mub");
  label_p->DrawLatex(0.2,0.87,"PYTHIA+HYDJET");

  profPanel_p->cd(4);

/*
  label_p->DrawLatex(.1, .92, Form("%s Lead Jet p_{T} > 120 GeV/c", algType[setNum]));
  label_p->DrawLatex(.1, .88, Form("%s Sublead Jet p_{T} > 50 GeV/c", algType[setNum]));
  label_p->DrawLatex(.1, .84, Form("%s Jet #Delta #phi > 5#pi/6", algType[setNum]));
  label_p->DrawLatex(.1, .80, Form("%s Lead/Sublead Jet |#eta| < .5", algType[setNum]));
  label_p->DrawLatex(.1, .76, Form("CMS Preliminary"));
*/
  profPanel_p->cd(3);
 
  makeHistForPtStack(getGraph2_p[0], getGraph2_p[1], getGraph2_p[2], getGraph2_p[3], getGraph2_p[4], getGraph2_p[5], hist2_p[0], hist2_p[1], hist2_p[2], hist2_p[3], hist2_p[4], hist2_p[5], 3);

  drawHistToPTStack(hist2_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(hist2_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[4], kRed + 1, "E1 HIST SAME");

  hist2_p[5]->DrawCopy("SAME E1");

  zeroLine_p->Draw();

  label_p->DrawLatex(.1,.87,"PYTHIA+HYDJET");
  label_p->DrawLatex(.1, .8, "0-30%");

  profPanel_p->cd(5);
  drawHistToPTStack(hist1_sub_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(hist1_sub_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(hist1_sub_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(hist1_sub_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(hist1_sub_p[4], kRed + 1, "E1 HIST SAME");

  hist1_sub_p[5]->DrawCopy("SAME E1");
  hist1_sub_p[5]->Write();
  
  zeroLine_p->Draw();
 
  label_p->DrawLatex(.2, .8, "30-100%");
  label_p->DrawLatex(.2,.87,"(PYTHIA+HYDJET)-PYTHIA");

  profPanel_p->cd(6);

  drawHistToPTStack(hist2_sub_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(hist2_sub_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(hist2_sub_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(hist2_sub_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(hist2_sub_p[4], kRed + 1, "E1 HIST SAME");

  hist2_sub_p[5]->DrawCopy("SAME E1");
  hist2_sub_p[5]->Write();
 
  zeroLine_p->Draw();

  label_p->DrawLatex(.2,.87,"(PYTHIA+HYDJET)-PYTHIA");
  label_p->DrawLatex(.2, .8, "0-30%");

  profPanel_p->Write();

  claverCanvasSaving(profPanel_p, Form("../pngDir/%s%sImbDelR%s%sPTStack_%s_%s", gorr, algType[setNum], perpProj, Corr, GLN, fileTag1), "png");

  delete label_p;
  delete zeroLine_p;
  delete profPanel_p;
  delete dummy_p;

  for(Int_t histIter = 0; histIter < 6; histIter++){
    delete hist1_p[histIter];
    delete hist2_p[histIter];
    delete hist_pp_p[histIter];
    delete hist1_sub_p[histIter];
    delete hist2_sub_p[histIter];
  }

  panelFile_p->Close();
  delete panelFile_p;
}

void cfmDiJet_DelRPlots_int(const char* inName,const char* inName_pp, const char* outName, Bool_t montecarlo = false)
{
  TH1::SetDefaultSumw2();

 if(!strcmp(inName, PbPbUp)){
    std::cout << PbPbUp << std::endl;
    fileTag1 = "PbPbUp";
  }
  else if(!strcmp(inName, PbPbDown)){
    std::cout << PbPbDown << std::endl;
    fileTag1 = "PbPbDown";
  }
  else if(!strcmp(inName, finalHydjet)){
    std::cout << finalHydjet << std::endl;
    fileTag1 = "finalHydjet";
  }
  else if(!strcmp(inName, PbPbData)){
    std::cout << PbPbData << std::endl;
    fileTag1 = "PbPbData";
  }


if(!strcmp(inName_pp, DataPP)){
    std::cout << DataPP << std::endl;
    fileTag2 = "DataPP";
  }
 else if(!strcmp(inName_pp, Pythia80)){
    std::cout << Pythia80 << std::endl;
    fileTag2 = "Pythia80";
  }
 else if(!strcmp(inName_pp, PPData_UP)){
    std::cout << PPData_UP << std::endl;
    fileTag2 = "PPData_UP";
  }
 else if(!strcmp(inName_pp, PPData_DOWN)){
    std::cout << PPData_DOWN << std::endl;
    fileTag2 = "PPData_DOWN";
  }
 else if(!strcmp(inName_pp, Pythia80_QG)){
    std::cout << Pythia80_QG << std::endl;
    fileTag2 = "Pythia80_QG";
  }
  else if(!strcmp(inName_pp, PythiaNew)){
    std::cout << PythiaNew << std::endl;
    fileTag2 = "PythiaNew";
  }
  else if(!strcmp(inName_pp, Pythia80_redo)){
    std::cout << Pythia80_redo << std::endl;
    fileTag2 = "Pythia80_redo";
  }
  else if(!strcmp(inName_pp, Pythia80_redoPF)){
    std::cout << Pythia80_redoPF << std::endl;
    fileTag2 = "Pythia80_redoPF";
  }
  else if(!strcmp(inName_pp, ppData_Down)){
    std::cout << ppData_Down << std::endl;
    fileTag2 = "ppData_Down";
  }

  std::cout << "Filetag1 is: " << fileTag1 << std::endl;
  std::cout << "Filetag2 is: " << fileTag2 << std::endl;

  inFile1_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile1_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  inFilepp_p = new TFile(inName_pp, "READ");
  inTreepp_p = (TTree*)inFilepp_p->Get("jetTree");
  inTreepp_p->AddFriend("trackTree");


  Int_t jetAlgMax = 4;

  if(montecarlo){
    inTree_p->AddFriend("genTree");
    jetAlgMax = 5;
    
    inTreepp_p->AddFriend("genTree");
  }

  const char* corr[1]={"Corr"};
  const char* ptBins[5] = {"0_1", "1_2", "2_4", "4_8", "8_100"};

  Int_t centLow[2] = {0, 60};
  Int_t centHi[2] = {59, 199};

  //These values are the Aj values x 100
  Float_t Aj_low[3] = {0,0,22};
  Float_t Aj_high[3] = {100,22,100}; 

  for(int ajIter = 0; ajIter<3; ajIter++){
    for(Int_t algIter = 3; algIter < jetAlgMax; algIter++){
      if(algIter != 3 && algIter !=4) continue;

      for(Int_t corrIter = 0; corrIter < 1; corrIter++){
        for(Int_t centIter = 0; centIter < 2; centIter++){
	  makeImbDelRGraph(inTree_p, outName, "r", algIter, "ProjA", "F", centLow[centIter], centHi[centIter], Aj_low[ajIter], Aj_high[ajIter], -10, 10, "N", corr[corrIter]);
          makeImbDelRGraphPP(inTreepp_p, outName, "r", algIter, "ProjA", "F",Aj_low[ajIter], Aj_high[ajIter], -10, 10, "N", corr[corrIter]);

	  if(montecarlo){      
	    makeImbDelRGraph(inTree_p, outName, "g", algIter, "ProjA", "F", centLow[centIter], centHi[centIter],Aj_low[ajIter], Aj_high[ajIter], -10, 10, "N", "");
            makeImbDelRGraphPP(inTreepp_p, outName, "g", algIter, "ProjA", "F",Aj_low[ajIter], Aj_high[ajIter], -10, 10, "N", "");
          }  

	  for(Int_t ptBinIter = 0; ptBinIter < 5; ptBinIter++){
	    makeImbDelRGraph(inTree_p, outName, "r", algIter, "ProjA", ptBins[ptBinIter], centLow[centIter], centHi[centIter],Aj_low[ajIter], Aj_high[ajIter], -10, 10, "N", corr[corrIter]);
            makeImbDelRGraphPP(inTreepp_p, outName, "r", algIter, "ProjA", ptBins[ptBinIter],Aj_low[ajIter], Aj_high[ajIter], -10, 10, "N", corr[corrIter]);
          
	    if(montecarlo){
	      makeImbDelRGraph(inTree_p, outName, "g", algIter, "ProjA", ptBins[ptBinIter], centLow[centIter], centHi[centIter],Aj_low[ajIter], Aj_high[ajIter], -10, 10, "N", "");
              makeImbDelRGraphPP(inTreepp_p, outName, "g", algIter, "ProjA", ptBins[ptBinIter],Aj_low[ajIter], Aj_high[ajIter], -10, 10, "N", "");
	      }
          
          }
        }
      }
    makeImbDelRPtStack(outName, "r", algIter,Aj_low[ajIter], Aj_high[ajIter], "ProjA", "N", "Corr");

    if(montecarlo)
      makeImbDelRPtStack(outName, "g", algIter,Aj_low[ajIter], Aj_high[ajIter], "ProjA", "N", "");
    }
  }


 inFile1_p->Close();
 delete inFile1_p;
 return;
}

