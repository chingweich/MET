#include <string>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TH1.h>
#include "TAttFill.h"
#include "TLine.h"
#include "TAxis.h"
using namespace std;
TFile*  f[7];
TH1D* h[9];
TCanvas* c1;
int COLOR[6]={1,14,42,44,46,48};
double factor[8]={1,1,19712.225*63.5/11764538,19712.225*39.4/12511326,19712.225*25.8/10783509,19712.225*56.0/9959752,19712.225*22.4/9910267,19712.225*7.6/9769891};
std::string dataset[6]={"data","uncertainties","DY_PTZ-70to100","DY_PTZ-100","Top","EWK"};
void drawSameBase(string ob,int aDraw){
  
  //Tag----------------------------------------------------------------
  double x1NDC = 0.7322;
  double y1NDC = 0.7835;
  double x2NDC = 0.9822;
  double y2NDC = 0.9860;
  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  leg->SetFillColor(18);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(1);
 //Tag----------------------------------------------------------------

//0fordata-----1foruncertainty-------2-7forMC-----------------------------------
  h[0]=(TH1D*) f[0]->FindObjectAny(ob.data());
  h[2]=(TH1D*) f[1]->FindObjectAny(ob.data());
  h[3]=(TH1D*) f[2]->FindObjectAny(ob.data());
  h[4]=(TH1D*) f[3]->FindObjectAny(ob.data());
  h[5]=(TH1D*) f[4]->FindObjectAny(ob.data());
  h[6]=(TH1D*) f[5]->FindObjectAny(ob.data());
  h[7]=(TH1D*) f[6]->FindObjectAny(ob.data());
  h[1]=(TH1D*) h[2]->Clone();
  
  //calculate uncertainty
  Int_t nBins=h[1]->GetNbinsX();
  double errH[6],errTotal;
  for(int i=0;i<nBins;i++){
	for(int j=0;j<6;j++){
		errH[j]=  h[j+2]->GetBinError(i);
		errH[j]*=factor[j+2];
        errTotal=errTotal+errH[j]*errH[j];
	}
	errTotal=sqrt(errTotal);
    h[1]->SetBinContent(i,errTotal);
  }  
   //setLogy
  c1->SetLogy();
  
  
  
  for(int i=0; i< 6; i++){
    if (i==5)cout<<ob.data()<<endl;
    //if(i==1)continue;
    //h[i] = (TH1F*)f->FindObjectAny(tempName[i].data());
    //h[i]->SetXTitle(xtitle.data());
    h[i]->Scale(factor[i]);
    h[i]->SetLineColor(COLOR[i]);
    h[i]->SetMarkerColor(COLOR[i]);
    h[i]->SetMarkerSize(1);
    //h[i]->SetMarkerStyle(MARKER[i]);
    //h[i]->SetMinimum(0);
    if(i==0){
      //c1->SetLogy();
      h[i]->Draw("e");
	  //c1->SetLogy();
    }
    else{
		//add Histo
		for(int j =i+1;j<8;j++){
			h[i]->Add(h[j],factor[j]);
		}
      
      //if(i==1)h[i]->Add(h[2],19712.225*39.4/12511326);
      h[i]->SetFillColor(COLOR[i]);
      h[i]->SetFillStyle(1001);
      h[i]->Draw("csame");
    }
    leg->AddEntry(h[i],dataset[i].data());
   
  }
  leg->Draw("same");
  h[0]->Draw("esame");
  cout<<ob.data()<<endl;

  
  if (aDraw==0)c1->Print("fig/met.pdf(");
  //else if (aDraw==4)c1->Print("fig/met.pdf)");
  else c1->Print("fig/met.pdf");
  c1->Print(Form("fig/%s.eps",ob.data()));

  c1->SetLogy(0);
  
  h[8]=(TH1D*) h[0]->Clone();
  h[1]->Add(h[2],-1);
  
  
  h[0]->Divide(h[2]);
  h[0]->SetYTitle("data/MC");
  double errData,errMc,averData,averMc,errRatio;
  for(int i=0;i<nBins;i++){
      errMc=h[1]->GetBinContent(i);
      errData=h[8]->GetBinError(i);
      averMc=h[2]->GetBinContent(i);
      averData=h[8]->GetBinContent(i);
      //cout<<i<<"=("<<errMc<<","<<errData<<","<<averMc<<","<<averData<<")"<<endl;
      errMc/=averMc;
      errData/=averData;
      if((averMc==0)&&(averData==0)){
          h[0]->SetBinError(i,0);
          continue;
      }
      //cout<<"("<<errMc<<","<<errData<<","<<averMc<<","<<averData<<")"<<endl;
      errRatio=sqrt(errMc*errMc+errData*errData);
      //cout<<"errRatio="<<errRatio<<endl;
      errRatio*=h[0]->GetBinContent(i);
      h[0]->SetBinError(i,errRatio);
  }
  

  h[0]->Draw("e");
  /*
  double systematic;
  for(int i=0;i<nBins;i++){
    systematic=(h[2]->GetBinContent(i)+h[1]->GetBinContent(i))/h[2]->GetBinContent(i);
    h[1]->SetBinContent(i,systematic);
  }
  h[1]->Draw("csame");
  for(int i=0;i<nBins;i++){
    systematic=(h[2]->GetBinContent(i)-h[1]->GetBinContent(i))/h[2]->GetBinContent(i);
    h[1]->SetBinContent(i,systematic);
  }
  h[1]->SetFillColor(0);
  h[1]->Draw("csame");
  h[0]->Draw("esame");
  */
  TAxis *xAxis = h[0]-> GetXaxis();
  double axisXMin=xAxis->GetXmin(),axisXMax=xAxis->GetXmax();
  
  TLine *ratio1 =new TLine(axisXMin,1,axisXMax,1);
  ratio1-> SetLineColor(2);
  ratio1->Draw("same");
  

  if(aDraw==4)c1->Print("fig/met.pdf)");
  else c1->Print("fig/met.pdf");

  
}


void drawSame(std::string s1,std::string s2,std::string s3,std::string s4,std::string s5,std::string s6,std::string s7){
   
  f[0]=new TFile(s1.data());
  f[1]=new TFile(s2.data());
  f[2]=new TFile(s3.data());
  f[3]=new TFile(s4.data());
  f[4]=new TFile(s5.data());
  f[5]=new TFile(s6.data());
  f[6]=new TFile(s7.data()); 
  string histoName[5]={"qT","uT","all","uParaAddqT","uPerp"};
  //TFile f2(s2.data());  
  //TFile f3(s3.data());
  //TH1D* qt[3],*ut[3],*all[3],*up[3],upp[3];
  
  
  std::string filename = "drawSame.pdf";
  

  gSystem->mkdir("fig");
  
  c1 = new TCanvas("c1","",700,500);
  gStyle->SetOptStat(0);

  c1->cd();
  
  for(int i=0 ;i<5 ;i++){
    drawSameBase(histoName[i],i);
    //c1->Print(Form("%s.pdf",histoName[i].data()));
  }
 
  //c1->Print(Form("%s(",filename.data()));

  //c1->Print(Form("%s)",filename.data()));

}














