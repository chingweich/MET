#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TH2D.h>
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
#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "setNCUStyle.C"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"


TCanvas* c1;

void MetScanning(){	
setNCUStyle();
	c1 = new TCanvas("c1","",800,600);
	TH1F* th1[10];
	th1[0]=new TH1F("pfMET","",20,0,2000);
	th1[1]=new TH1F("pfMETCut","",20,0,2000);
	th1[2]=new TH1F("caloMET","",100,0,2000);
	th1[3]=new TH1F("caloMETCut","",100,0,2000);
	th1[4]=new TH1F("pfClusterMETPt","",100,0,2000);
	th1[5]=new TH1F("pfClusterMETPtCut","",100,0,2000);
	th1[6]=new TH1F("pfCaloMETPt","",100,0,2000);
	th1[7]=new TH1F("pfCaloMETPtCut","",100,0,2000);
	th1[0]->SetXTitle("pfMET[GeV]");
	th1[2]->SetXTitle("caloMET[GeV]");
	th1[4]->SetXTitle("pfClusterMET[GeV]");
	th1[6]->SetXTitle("pfCaloMET[GeV]");
	
	for(int i=0;i<8;i++)th1[i]->SetYTitle("Events");
	
	TH2D* th2[10];
	th2[0]=new TH2D("caloMET_vs_pfMET","",100,0,2000,100,0,2000);
	th2[0]->SetXTitle("pfMET[GeV]");th2[0]->SetYTitle("caloMET[GeV]");
	
	th2[1]=new TH2D("pfClusterMET_vs_pfMET","",100,0,2000,100,0,2000);
	th2[1]->SetXTitle("pfMET[GeV]");th2[1]->SetYTitle("pfClusterMET[GeV]");
	
	th2[2]=new TH2D("pfCaloMET_vs_pfMET","",100,0,2000,100,0,2000);
	th2[2]->SetXTitle("pfMET[GeV]");th2[2]->SetYTitle("pfCaloMET[GeV]");
	
	th2[3]=new TH2D("caloMET_vs_pfMETC","",100,0,2000,100,0,2000);
	th2[3]->SetXTitle("pfMET[GeV]");th2[3]->SetYTitle("caloMET[GeV]");
	
	th2[4]=new TH2D("pfClusterMET_vs_pfMETC","",100,0,2000,100,0,2000);
	th2[4]->SetXTitle("pfMET[GeV]");th2[4]->SetYTitle("pfClusterMET[GeV]");
	
	th2[5]=new TH2D("pfCaloMET_vs_pfMETC","",100,0,2000,100,0,2000);
	th2[5]->SetXTitle("pfMET[GeV]");th2[5]->SetYTitle("pfCaloMET[GeV]");
	
	//th2[0]->SetTitle("beforeFilter");th2[1]->SetTitle("beforeFilter");th2[2]->SetTitle("beforeFilter");
	//th2[3]->SetTitle("afterAllFilter");th2[4]->SetTitle("afterAllFilter");th2[5]->SetTitle("afterAllFilter");
	TFile *f;
	TTree *tree;
	
	
	string  st="/data7/chchen/METScanning/HotLineBCDE.root";
	//string  st="tuple_5.root";
	f=TFile::Open(Form("%s",st.data()));
	//TDirectory * dir;
	 //dir = (TDirectory*)f->Get(Form("%s:/",st.data()));
		
	f->GetObject("tree",tree);
	TreeReader data(tree);
	data.Print();
	cout<<"Entry="<<data.GetEntriesFast()<<endl;
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
		if((jEntry%1000)==0)cout<<"processing"<<1000*jEntry<<endl;
		data.GetEntry(jEntry);
		Float_t pfMETPt = data.GetFloat("pfMETPt");
		Float_t  pfClusterMETPt = data.GetFloat("pfClusterMETPt");
		Float_t  pfCaloMETPt = data.GetFloat("pfCaloMETPt");
		Float_t  caloMETPt = data.GetFloat("caloMETPt");
                Long64_t  run = data.GetLong64("run");
                if(run<274000)continue;
		
		//cout<<pfMETPt[0]<<endl;
		th1[0]->Fill(pfMETPt);
		th1[2]->Fill(caloMETPt);
		th1[4]->Fill(pfClusterMETPt);
		th1[6]->Fill(pfCaloMETPt);
		
		th2[0]->Fill(pfMETPt,caloMETPt);
		th2[1]->Fill(pfMETPt,pfClusterMETPt);
		th2[2]->Fill(pfMETPt,pfCaloMETPt);
		
		Bool_t   filter_csc2015=data.GetBool  ("filter_csc2015"); 
		Bool_t   filter_globaltighthalo2016=data.GetBool  ("filter_globaltighthalo2016"); 
		Bool_t   filter_globalsupertighthalo2016=data.GetBool  ("filter_globalsupertighthalo2016"); 
		Bool_t   filter_hcalstriphalo=data.GetBool  ("filter_hcalstriphalo"); 
		Bool_t   filter_hbher2l=data.GetBool  ("filter_hbher2l"); 
		Bool_t   filter_hbher2t=data.GetBool  ("filter_hbher2t"); 
		Bool_t   filter_hbheiso=data.GetBool  ("filter_hbheiso"); 
		Bool_t   filter_ecaltp=data.GetBool  ("filter_ecaltp"); 
		Bool_t   filter_ecalsc=data.GetBool  ("filter_ecalsc"); 
            Bool_t   filter_badChCand=data.GetBool  ("filter_badChCand"); 
		Bool_t   filter_badPFMuon=data.GetBool  ("filter_badPFMuon"); 
		//if(!filter_csc2015)continue;
		if(!filter_globaltighthalo2016)continue;

		if(!filter_hbher2l)continue;

		if(!filter_hbheiso)continue;

		if(!filter_ecaltp)continue;

		if(!filter_ecalsc)continue;

            if(!filter_badChCand)continue;

		if(!filter_badPFMuon)continue;
		
		th1[1]->Fill(pfMETPt);
		th1[3]->Fill(caloMETPt);
		th1[5]->Fill(pfClusterMETPt);
		th1[7]->Fill(pfCaloMETPt);
		
		th2[3]->Fill(pfMETPt,caloMETPt);
		th2[4]->Fill(pfMETPt,pfClusterMETPt);
		th2[5]->Fill(pfMETPt,pfCaloMETPt);
		
		
	}
		th1[1]->SetLineColor(2);
		th1[3]->SetLineColor(2);
		th1[5]->SetLineColor(2);
		th1[7]->SetLineColor(2);
		
		TLegend* leg ;
		leg=new TLegend(0.691452,0.762447,0.890645,0.883966);
		leg->AddEntry(th1[0],"beforeFilter");
		leg->AddEntry(th1[1],"afterAllFilter");
		 leg->SetTextSize(0.04);
		TLatex *lar = new TLatex();

		lar->SetNDC(kTRUE);
		lar->SetTextSize(0.07);
		lar->SetLineWidth(5);
		lar->SetTextAlign(14);
		//lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
		
		
		TLatex *lar2 = new TLatex();

		lar2->SetNDC(kTRUE);
		lar2->SetTextSize(0.04);
		lar2->SetLineWidth(5);
		lar2->SetTextAlign(12);
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		/* 
		c1->SetLogy(1);
		th1[0]->Draw("");
		th1[1]->Draw("same");
		leg->Draw("same");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/pfMETPt.pdf");
		
		th1[2]->Draw("");
		th1[3]->Draw("same");
		leg->Draw("same");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/caloMETPt.pdf");
		
		
		th1[4]->Draw("");
		th1[5]->Draw("same");
		leg->Draw("same");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/pfClusterMETPt.pdf");
		
		th1[6]->Draw("");
		th1[7]->Draw("same");
		leg->Draw("same");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/pfCaloMETPt.pdf");
		
		c1->SetLogy(0);
		th2[0]->Draw("colz");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/caloMET_vs_pfMET.pdf");
		th2[1]->Draw("colz");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/pfClusterMET_vs_pfMET.pdf");
		th2[2]->Draw("colz");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/pfCaloMET_vs_pfMET.pdf");
		
		th2[3]->Draw("colz");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/caloMET_vs_pfMETCut.pdf");
		th2[4]->Draw("colz");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/pfClusterMET_vs_pfMETCut.pdf");
		th2[5]->Draw("colz");
		lar->DrawLatex(0.14, 0.94, "CMS");
		lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		c1->Print("plot/pfCaloMET_vs_pfMETCut.pdf");
		*/
		
		th1[2]=new TH1F("pfMET","",21,0,2100);
		th1[3]=new TH1F("pfMETCut","",21,0,2100);
		
		for(int i=0;i<21;i++){
			th1[2]->SetBinContent(i+1,th1[0]->GetBinContent(i+1));
			th1[3]->SetBinContent(i+1,th1[1]->GetBinContent(i+1));
		}
		
		TGraphAsymmErrors * tga1= new TGraphAsymmErrors(th1[3],th1[2]);
		
		//th1[0]->SetBinContent(21,th1[0]->GetBinContent(21));
		//th1[1]->SetBinContent(21,th1[1]->GetBinContent(21));
		tga1->GetXaxis()->SetNdivisions(508);
		tga1->GetXaxis()->SetRangeUser(0,2100);
		tga1->GetXaxis()->SetTitle("E^{miss}_{T}[GeV]");
		tga1->GetYaxis()->SetTitle("tagging rate");
		tga1->SetLineColor(1);
		tga1->SetMarkerSize(tga1->GetMarkerSize()*1.5);
		tga1->Draw("AP");
		
		
		lar2->DrawLatex(0.66, 0.94, "15.9 fb^{-1} (13 TeV,2016)");
		//lar2->SetTextAlign(18);
		
		lar->SetTextAlign(12);
		lar->DrawLatex(0.7, 0.7, "CMS");
		lar2->SetTextAlign(11);
		lar2->DrawLatex(0.7, 0.63, "#it{#bf{Preliminary}}");
		c1->SetLogy(0);
		//lar->DrawLatex(0.14, 0.94, "CMS");
		//lar2->DrawLatex(0.25, 0.94, "Preliminary               15.9fb-1(13 TeV,2016)");
		
		c1->Print("plot/METRatio.pdf");
		
}
