#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>

TFile *f;
TTree *tree;

void drawMetFig(){
	gStyle->SetOptStat(0000000000);
		f = TFile::Open("METScanningHotLine2015C.root");	
		TDirectory * dir = (TDirectory*)f->Get("METScanningHotLine2015C.root:/Metanalyzer");
		
		TH1F* pfmetNofilter= (TH1F*) dir->Get("h_pfMETNoFilter");
		TH1F* pfmethbhet= (TH1F*) dir->Get("h_pfMEThbhet");
		TH1F* pfmetcsct= (TH1F*) dir->Get("h_pfMETcsct");
		TH1F* pfmeteesc= (TH1F*) dir->Get("h_pfMETEESCBad");
		TH1F* pfmet2015= (TH1F*) dir->Get("h_pfMETCSCTH2015");
		TH1F* pfmetmu= (TH1F*) dir->Get("h_pfMETCSCTHMu");
		TH1F* pfmethcal= (TH1F*) dir->Get("h_pfMETHCALStrp");
		TH1F* pfmetECDC= (TH1F*) dir->Get("h_pfMETECDCTP");
		TH1F* pfmet2015a= (TH1F*) dir->Get("h_pfMETAllCSCTH2015");
		TH1F* pfmetmua= (TH1F*) dir->Get("h_pfMETAllCSCTHMu");
		TH1F* pfmethcala= (TH1F*) dir->Get("h_pfMETAllHCALStrp");
		TH1F* pfmethcala2015= (TH1F*) dir->Get("h_pfMETAllCSCTH2015HCALStrp");
		
		pfmetNofilter->Rebin(10);
		pfmethbhet->Rebin(10);
		pfmetcsct->Rebin(10);
		pfmeteesc->Rebin(10);
		pfmet2015->Rebin(10);
		pfmetmu->Rebin(10);
		pfmethcal->Rebin(10);
		pfmetECDC->Rebin(10);
		pfmet2015a->Rebin(10);
		pfmetmua->Rebin(10);
		pfmethcala->Rebin(10);
		pfmethcala2015->Rebin(10);
		
		TCanvas* c1;
		
		c1 = new TCanvas("c1","",1360,768);
		c1->SetLogx(1);
		c1->SetLogy(1);
		double x1NDC = 0.6522;
        double y1NDC = 0.7835;
		double x2NDC = 0.9122;
		double y2NDC = 0.9860;
		TLegend* leg ,*leg2 ;
    leg=new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
	leg->SetFillColor(18);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(2);
    //leg->AddEntry(th23,"signal-600to4500");
	
	//leg->AddEntry(th24,"DY");
	leg->AddEntry(pfmetNofilter,"Nofilter");
	leg->AddEntry(pfmethbhet,"hbhet");
	leg->AddEntry(pfmetcsct,"csct");
	leg->AddEntry(pfmeteesc,"eebadscf");
	leg->AddEntry(pfmetECDC,"ecaldeadcellTPFilter");
	
	
	pfmetNofilter->SetLineColor(1);
	pfmethbhet->SetLineColor(2);
	pfmetcsct->SetLineColor(3);
	pfmeteesc->SetLineColor(4);
	pfmetECDC->SetLineColor(6);
	
	
	
	pfmetNofilter->Draw();
	pfmetcsct->Draw("same");
	pfmeteesc->Draw("same");
	pfmetECDC->Draw("same");
	pfmethbhet->Draw("same");	
	
	leg->Draw("same");
    c1->Print("pdf1.png");
	
	
	leg->Clear();
	
	leg->AddEntry(pfmet2015,"CSCTightHalo2015Filter");
	leg->AddEntry(pfmetmu,"CSCTightHaloTrkMuUnvetoFilter");
	leg->AddEntry(pfmethcal,"HcalStripHaloFilter");
	
	pfmethcal->SetLineColor(2);
	pfmetmu->SetLineColor(3);
	pfmet2015->SetLineColor(4);
	
	pfmethcal->Draw();
	pfmetmu->Draw("same");
	pfmet2015->Draw("same");
	leg->Draw("same");
	c1->Print("pdf2.png");
	
	
	leg->Clear();
	leg2=new TLegend(0.55,y1NDC,x2NDC,y2NDC);
	leg2->SetFillColor(18);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.025);
    leg2->SetBorderSize(2);
	leg2->AddEntry(pfmet2015a,"hbhet+csct+eebadscf+CSCTightHalo2015Filter");
	leg2->AddEntry(pfmetmua,"hbhet+csct+eebadscf+CSCTightHaloTrkMuUnvetoFilter");
	leg2->AddEntry(pfmethcala,"hbhet+csct+eebadscf+HcalStripHaloFilter");
	leg2->AddEntry(pfmethcala2015,"All + CSC2015 + HCAL strip");
	
	pfmethcala->SetLineColor(2);
	pfmetmua->SetLineColor(3);
	pfmet2015a->SetLineColor(4);
	pfmethcala2015->SetLineColor(6);
	
	
	pfmetmua->Draw();
	pfmethcala->Draw("same");
	pfmet2015a->Draw("same");
	pfmethcala2015->Draw("same");
	leg2->Draw("same");
	c1->Print("pdf4.png");
	
	
}