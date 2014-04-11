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
TTree* inTree_p = 0;

TFile* outFile_p = 0;

const char* algType[5] = {"PuPF", "PuCalo", "VsPF", "VsCalo", "T"};

//append to every histo so know which sample via shorthand on workblog                                                          

const char* fileTag1;

//shorthands, w/ _CFMSKIM.h                                                                                                                      

const char* Di80a = "Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMSKIM.root";

const char* Di80b = "Dijet80_HydjetDrum_v27_mergedV1_CFMSKIM.root";

const char* Di80c = "HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v2_CFMSKIM.root";

const char* Di80d = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";                            

const char* Di80e = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet14_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";

const char* Di80f = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_50k_CFMSKIM.root";

const char* Di80g = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";

const char* Di80h = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";

const char* Di100a = "Dijet100_HydjetDrum_v27_mergedV1_CFMSKIM.root";

const char* EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80_CFMSKIM.root";

const char* DataA = "Track8_Jet17_GR_R_53_LV6_SUB_0_CFMSKIM.root";

const char* DataB = "hiForest_Jet80or95_GR_R_53_LV6_02Mar2014_1300CET_Track8_Jet15_0_1200k_CFMSKIM.root";

const char* DataC = "hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21_0_700k_CFMSKIM.root";

const char* DataD = "hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0_1200k_CFMSKIM.root";

const char* DataE = "HydjetDrum_Pyquen_Dijet100_Embedded_FOREST_STARTHI53_LV1_Track8_Jet21_290k_v0_merged_0_CFMSKIM_20140322_0.root";

const char* DataF = "HydjetDrum_Pyquen_Dijet120_Embedded_FOREST_STARTHI53_LV1_Track8_Jet21_290k_v0_merged_0_CFMSKIM_20140322_0.root";

const char* DataG = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM_avgaxis.root";

const char* DataH = "hiForest_Jet80or95_GR_R_53_Lv6_12Mar2014_000CET_Track8_Jet21_0_1200k.root";

const char* DataI = "HydjetDrum_Pyquen_Dijet120_Embedded_FOREST_STARTHI53_LV1_Track8_Jet24_290k_v0_merged_0_CFMSKIM_20140322_0.root";

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


void makeImbDelRGraph(TTree* getTree_p, const char* outName, const char* gorr, Int_t setNum, const char* perpProj, const char* FPT, Int_t centLow, Int_t centHi, Int_t graphLow, Int_t graphHi, const char* GLN = "N", const char* Corr = "")
{
  inFile1_p->cd();

  Int_t setCorrNum = setNum;
  if(strcmp("", Corr) != 0)
    setCorrNum = setNum + 5;

  const char* title = Form("%s%sImbDelR%s%s%s_%d%d_%s_%s_g", gorr, algType[setNum], perpProj, FPT, Corr, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), GLN, fileTag1);

  //const char* title_integrated = Form("%s%sImbDelR%s%s%s_%d%d_%s_%s_integrated_g", gorr, algType[setNum], perpProj, FPT, Corr, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), GLN, fileTag1);

  TGraphErrors* imbDelRGraph_p = new TGraphErrors(10);
  imbDelRGraph_p->GetXaxis()->SetLimits(0.00, 2.00);
  niceTGraphErrors(imbDelRGraph_p, graphHi, graphLow);

  //TGraphErrors* netimbDelRGraph_p = new TGraphErrors(10);
  //netimbDelRGraph_p->GetXaxis()->SetLimits(0.00, 2.00);
  //niceTGraphErrors(netimbDelRGraph_p, graphHi, graphLow);

  TH1F* getHist_p;

  const char* coneR[10] = {"1C", "2C", "3C", "4C", "5C", "6C", "7C", "8C", "9C", "10C"};

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 0.5, GLN);
  TCut phiCut = makeDelPhiCut(setNum,5*TMath::Pi()/6.0);
  TCut assymCut = makeAsymmCut(setNum,0.0,1);
  TCut jetLCut = Form("AlgLeadJtPt[%d] > 120", setNum);

  const char* name1[10] = {"1C(10000, -10000, 10000)", "2C(10000, -10000, 10000)", "3C(10000, -10000, 10000)", "4C(10000, -10000, 10000)", "5C(10000, -10000, 10000)", "6C(10000, -10000, 10000)", "7C(10000, -10000, 10000)", "8C(10000, -10000, 10000)", "9C(10000, -10000, 10000)", "10C(10000, -10000, 10000)"};

  Float_t point[10] = {.1, .3, .5, .7, .9, 1.1, 1.3, 1.5, 1.7, 1.9};
  Float_t xErr[10] = {.1, .1, .1, .1, .1, .1, .1, .1, .1, .1};
  Float_t netimb=0.0; 
  Float_t netimberr=0.0; 

  for(Int_t binIter = 0; binIter < 10; binIter++){
    TString var = Form("%sAlgImb%s%s%s[%d]", gorr, perpProj, coneR[binIter], FPT, setCorrNum);
    getTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && assymCut);

    getHist_p = (TH1F*)inFile1_p->Get(coneR[binIter]);
    netimb += getHist_p->GetMean();
    netimberr += TMath::Power(getHist_p->GetMeanError(),2);
    //netimberr += getHist_p->GetMeanError();

    //imbDelRGraph_p->SetPoint(binIter, point[binIter], getHist_p->GetMean());
    //imbDelRGraph_p->SetPointError(binIter, xErr[binIter], getHist_p->GetMeanError());

    imbDelRGraph_p->SetPoint(binIter, point[binIter], netimb);
    imbDelRGraph_p->SetPointError(binIter, xErr[binIter], TMath::Power(netimberr,0.5));
    //imbDelRGraph_p->SetPointError(binIter, xErr[binIter], netimberr);

    //netimbDelRGraph_p->SetPoint(binIter, point[binIter], netimb);
    //netimbDelRGraph_p->SetPoint(binIter, point[binIter], netimberr);
    
  }

  imbDelRGraph_p->GetXaxis()->SetLimits(0.00, 2.00);
  niceTGraphErrors(imbDelRGraph_p, graphHi, graphLow);

 // netimbDelRGraph_p->GetXaxis()->SetLimits(0.00, 2.00);
  //niceTGraphErrors(netimbDelRGraph_p, graphHi, graphLow);

  outFile_p = new TFile(outName, "UPDATE");
  imbDelRGraph_p->Write(title);
  //netimbDelRGraph_p->Write(title_integrated);
  outFile_p->Close();

  delete outFile_p;
  delete imbDelRGraph_p;
  //delete netimbDelRGraph_p;
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

    h0_1_p->SetBinContent(iter + 1, sumYForPTStack(y0_1[iter], y1_2[iter], y2_4[iter], y4_8[iter], y8_100[iter]));
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



void drawHistToPTStack(TH1F* drawHist_p, Int_t color, const char* drawOpt)
{
  drawHist_p->SetFillColor(color);
  drawHist_p->SetMarkerStyle(6);
  drawHist_p->SetMarkerSize(.5);
  drawHist_p->Draw(drawOpt);
  drawHist_p->Draw("E1 SAME");
}



void makeImbDelRPtStack(const char* fileName, const char* gorr, Int_t setNum, const char* perpProj, const char* GLN, const char* Corr = "")
{
  TFile* panelFile_p = new TFile(fileName, "UPDATE");

  TGraphErrors* getGraph1_p[6];
  TGraphErrors* getGraph2_p[6];

  TH1F* hist1_p[6];
  TH1F* hist2_p[6];

  Float_t binArrayX[11] = {.00, .20, .40, .60, .80, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

  const char* namePT[6] = {"0_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 0; histIter < 6; histIter++){
    hist1_p[histIter] = new TH1F(Form("hist1%s_p", namePT[histIter]), Form("hist1%s_p", namePT[histIter]), 10, binArrayX);
    hist2_p[histIter] = new TH1F(Form("hist2%s_p", namePT[histIter]), Form("hist2%s_p", namePT[histIter]), 10, binArrayX);

    niceTH1(hist1_p[histIter], 60, -60, 505, 406);
    niceTH1(hist2_p[histIter], 60, -60, 505, 406);

    getGraph1_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbDelR%s%s%s_30100_%s_%s_g", gorr, algType[setNum], perpProj, namePT[histIter], Corr, GLN, fileTag1));
    getGraph2_p[histIter] = (TGraphErrors*)panelFile_p->Get(Form("%s%sImbDelR%s%s%s_030_%s_%s_g", gorr, algType[setNum], perpProj, namePT[histIter], Corr, GLN, fileTag1));
  }

  makeHistForPtStack(getGraph1_p[0], getGraph1_p[1], getGraph1_p[2], getGraph1_p[3], getGraph1_p[4], getGraph1_p[5], hist1_p[0], hist1_p[1], hist1_p[2], hist1_p[3], hist1_p[4], hist1_p[5], 1);

  TCanvas* profPanel_p;
  profPanel_p = new TCanvas(Form("%s%sImbDelR%s%sPTStack_%s_%s_c", gorr, algType[setNum], perpProj, Corr, GLN, fileTag1), Form("%s%sImbDelR%s%sPTStack_%s_%s_c", gorr, algType[setNum], perpProj, Corr, GLN, fileTag1), 700, 500);
  profPanel_p->Divide(2, 1, 0, 0);

  TLegend* leg;
  if(strcmp(gorr, "g") == 0)
    //leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("Truth #slash{p}_{T}^{||} v. %s #Delta R, %s, A_{J}>0.22", algType[setNum], fileTag1));
  //else
    //leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("%sTrk #slash{p}_{T}^{||} v. %s #Delta R, %s, A_{J}>0.22", Corr, algType[setNum], fileTag1));

leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("Truth #slash{p}_{T}^{||} v. %s #Delta R, %s", algType[setNum], fileTag1));
  else
    leg = new TLegend(0.15, 0.75, 0.95, 0.95, Form("%sTrk #slash{p}_{T}^{||} v. %s #Delta R, %s", Corr, algType[setNum], fileTag1));


  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(.04);
  leg->SetBorderSize(0);

  profPanel_p->cd(1);

  drawHistToPTStack(hist1_p[0], kBlue - 9, "E1 HIST");
  leg->AddEntry(hist1_p[0], ".5 < p_{T} < 1", "F");
  drawHistToPTStack(hist1_p[1], kYellow - 9, "E1 HIST SAME");
  leg->AddEntry(hist1_p[1], "1 < p_{T} < 2", "F");
  drawHistToPTStack(hist1_p[2], kOrange + 1, "E1 HIST SAME");
  leg->AddEntry(hist1_p[2], "2 < p_{T} < 4", "F");
  drawHistToPTStack(hist1_p[3], kGreen + 3, "E1 HIST SAME");
  leg->AddEntry(hist1_p[3], "4 < p_{T} < 8", "F");
  drawHistToPTStack(hist1_p[4], kRed + 1, "E1 HIST SAME");
  leg->AddEntry(hist1_p[4], "8 < p_{T}", "F");

  hist1_p[5]->Draw("SAME E1");
  leg->Draw("SAME");

  TLine* zeroLine_p = new TLine(0., 0., 2.0, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->DrawLatex(.2, .3, "30-100%");
  
  profPanel_p->cd(2);

  makeHistForPtStack(getGraph2_p[0], getGraph2_p[1], getGraph2_p[2], getGraph2_p[3], getGraph2_p[4], getGraph2_p[5], hist2_p[0], hist2_p[1], hist2_p[2], hist2_p[3], hist2_p[4], hist2_p[5], 2);

  drawHistToPTStack(hist2_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(hist2_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(hist2_p[4], kRed + 1, "E1 HIST SAME");

  hist2_p[5]->Draw("SAME E1");

  zeroLine_p->Draw();

  label_p->DrawLatex(.2, .3, "0-30%");

  label_p->DrawLatex(.1, .92, Form("%s Lead Jet p_{T} > 120 GeV/c", algType[setNum]));
  label_p->DrawLatex(.1, .88, Form("%s Sublead Jet p_{T} > 50 GeV/c", algType[setNum]));
  label_p->DrawLatex(.1, .84, Form("%s Jet #Delta #phi > 5#pi/6", algType[setNum]));
  label_p->DrawLatex(.1, .80, Form("%s Lead/Sublead Jet abs(#eta) < .5", algType[setNum]));

  profPanel_p->Write();

  claverCanvasSaving(profPanel_p, Form("../pngDir/%s%sImbDelR%s%sPTStack_%s_%s", gorr, algType[setNum], perpProj, Corr, GLN, fileTag1), "png");

  delete label_p;
  delete zeroLine_p;
  delete profPanel_p;

  for(Int_t histIter = 0; histIter < 6; histIter++){
    delete hist1_p[histIter];
    delete hist2_p[histIter];
  }

  panelFile_p->Close();
  delete panelFile_p;
}

void cfmDiJet_DelRPlots_int(const char* inName, const char* outName, Bool_t montecarlo = false)
{
  TH1::SetDefaultSumw2();

  if(!strcmp(inName, Di80a)){
    std::cout << Di80a << std::endl;
    fileTag1 = "Di80a";
  }
  else if(!strcmp(inName, Di80b)){
    std::cout << Di80b << std::endl;
    fileTag1 = "Di80b";
  }
  else if(!strcmp(inName, Di100a)){
    std::cout << Di100a << std::endl;
    fileTag1 = "Di100a";
  }
  else if(!strcmp(inName, EmDi80a)){
    std::cout << EmDi80a << std::endl;
    fileTag1 = "EmDi80a";
  }
  else if(!strcmp(inName, Di80c)){
    std::cout << Di80c << std::endl;
    fileTag1 = "Di80c";
  }
  else if(!strcmp(inName, Di80d)){
    std::cout << Di80d << std::endl;
    fileTag1 = "Di80d";
  }
  else if(!strcmp(inName, Di80e)){
    std::cout << Di80e << std::endl;
    fileTag1 = "Di80e";
  }
  else if(!strcmp(inName, Di80f)){
    std::cout << Di80f << std::endl;
    fileTag1 = "Di80f";
  }
  else if(!strcmp(inName, Di80g)){
    std::cout << Di80g << std::endl;
    fileTag1 = "Di80g";
  }
  else if(!strcmp(inName, Di80h)){
    std::cout << Di80h << std::endl;
    fileTag1 = "Di80h";
  }
  else if(!strcmp(inName, DataA)){
    std::cout << DataA << std::endl;
    fileTag1 = "DataA";
  }
  else if(!strcmp(inName, DataB)){
    std::cout << DataB << std::endl;
    fileTag1 = "DataB";
  }
  else if(!strcmp(inName, DataC)){
    std::cout << DataC << std::endl;
    fileTag1 = "DataC";
  }
  else if(!strcmp(inName, DataD)){
    std::cout << DataD << std::endl;
    fileTag1 = "DataD";
  }
  else if(!strcmp(inName, DataE)){
    std::cout << DataE << std::endl;
    fileTag1 = "DataE";
  }
  else if(!strcmp(inName, DataF)){
    std::cout << DataF << std::endl;
    fileTag1 = "DataF";
  }
   else if(!strcmp(inName, DataG)){
    std::cout << DataG << std::endl;
    fileTag1 = "DataG";
  }
   else if(!strcmp(inName, DataH)){
    std::cout << DataH << std::endl;
    fileTag1 = "DataH";
  }
 else if(!strcmp(inName, DataI)){
    std::cout << DataI << std::endl;
    fileTag1 = "DataI";
  }


  std::cout << "Filetag1 is: " << fileTag1 << std::endl;

  inFile1_p = new TFile(inName, "READ");
  inTree_p = (TTree*)inFile1_p->Get("jetTree");
  inTree_p->AddFriend("trackTree");

  Int_t jetAlgMax = 4;

  if(montecarlo){
    inTree_p->AddFriend("genTree");
    jetAlgMax = 5;
  }

  //const char* corr[2] = {"", "Corr"};
  const char* corr[1]={"Corr"};
  const char* ptBins[5] = {"0_1", "1_2", "2_4", "4_8", "8_100"};

  //  Int_t centLow[4] = {0, 20, 60, 100};
  //  Int_t centHi[4] = {19, 59, 99, 199};

  Int_t centLow[2] = {0, 60};
  Int_t centHi[2] = {59, 199};



  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){
    if(algIter != 3 && algIter !=4) continue;

    for(Int_t corrIter = 0; corrIter < 1; corrIter++){
      for(Int_t centIter = 0; centIter < 2; centIter++){
	makeImbDelRGraph(inTree_p, outName, "r", algIter, "ProjA", "F", centLow[centIter], centHi[centIter], -10, 10, "N", corr[corrIter]);
        makeImbDelRGraph(inTree_p, outName, "r", algIter, "ProjA", "F", centLow[centIter], centHi[centIter], -10, 10, "N", "");

	if(montecarlo)      
	  makeImbDelRGraph(inTree_p, outName, "g", algIter, "ProjA", "F", centLow[centIter], centHi[centIter], -10, 10, "N", "");

	for(Int_t ptBinIter = 0; ptBinIter < 5; ptBinIter++){
	  makeImbDelRGraph(inTree_p, outName, "r", algIter, "ProjA", ptBins[ptBinIter], centLow[centIter], centHi[centIter], -10, 10, "N", corr[corrIter]);
          makeImbDelRGraph(inTree_p, outName, "r", algIter, "ProjA", ptBins[ptBinIter], centLow[centIter], centHi[centIter], -10, 10, "N", "");

	  if(montecarlo)
	    makeImbDelRGraph(inTree_p, outName, "g", algIter, "ProjA", ptBins[ptBinIter], centLow[centIter], centHi[centIter], -10, 10, "N", "");

	}
      }
    }

    // makeImbDelRPtStack(outName, "r", algIter, "ProjA", "N");
    makeImbDelRPtStack(outName, "r", algIter, "ProjA", "N", "Corr");
    makeImbDelRPtStack(outName, "r", algIter, "ProjA", "N", "");

    if(montecarlo)
      makeImbDelRPtStack(outName, "g", algIter, "ProjA", "N", "");
    
  }



 inFile1_p->Close();
 delete inFile1_p;
 return;
}

