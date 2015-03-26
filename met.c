#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <TH2.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <cmath>
#include "TCanvas.h"
#include "TLine.h"
#include "untuplizer.h"
#include "TMarker.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"




void met(std::string fin){
  std::vector<string> infiles;
  //TString outputFileName;
  bool readOneFile=true;
  if(fin.find(".txt")!= std::string::npos)
    {
      readOneFile=false;
      TString endfix=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"${test%%.txt*}\"",fin.data()));
      //outputFileName = Form("histo_oot_%s_%d_timeW%02i.root",endfix.Data(),processType,(int)readoutWindow);
      //cout << "Output file = " << outputFileName << endl;
    }
  
  if(readOneFile){
      infiles.push_back(fin);
      TString endfix=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"${test%%.root*}\"",fin.data()));
      //outputFileName = Form("absHisto_oot_%s_%d_timeW%02i.root",endfix.Data(),processType,(int)readoutWindow);
      //cout << "Output file = " << outputFileName << endl;
    }
  
  else{
      FILE *fTable = fopen(fin.data(),"r");
      int flag=1;
      int nfile=0;

      while (flag!=-1){
         char filename[300];
         flag=fscanf(fTable,"%s",filename);
         // first reading input file
         if(flag!=-1){
       	 std::string tempFile = filename;    
     	 infiles.push_back(tempFile);
 	 nfile++;
         }
      }
     cout << "nfiles = " << nfile << endl;
  }

  TH1D* qT = new TH1D("qT", "Removed Vector Boson", 200, 0, 200);
  TH1D* uT = new TH1D("uT", "Sum of All Particels except Removed Vector Boson ", 100, 0, 200);
  TH1D* all = new TH1D("all", "Sum of All Particels ", 400, 0, 200);
  TH1D* uP = new  TH1D(" uParaAddqT", "Additon of uT parallel and qT", 600, -300, 300);
  TH1D* uPp = new  TH1D(" uPerp", "uT perpendicular ", 600, -300, 300);


  TreeReader data(infiles); // v5.3.12
  data.Print();
  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {


     // for (Long64_t ev = 0; ev < 100; ev++) {
     // print progress
     if (ev % 90000 == 0){
        fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());
     }
     data.GetEntry(ev);
   
     TLorentzVector u;
     u.SetPtEtaPhiM(0,0,0,0);
     TLorentzVector q;
     q.SetPtEtaPhiM(0,0,0,0);

     bool flag=0;

     Int_t nJet= data.GetInt("AK5nJet"); 
     Int_t nMu= data.GetInt("nMu"); 
     Int_t nEl= data.GetInt("nEle");
     Int_t nGen =data.GetInt("nGenPar");
     Int_t*  genId =data.GetPtrInt("genParId");
     Float_t* elePt    = data.GetPtrFloat("elePt");
     Float_t* eleEta   = data.GetPtrFloat("eleEta");
     Float_t* elePhi   = data.GetPtrFloat("elePhi");
     Float_t* muPt    = data.GetPtrFloat("muPt");
     Float_t* muEta   = data.GetPtrFloat("muEta");
     Float_t* muPhi   = data.GetPtrFloat("muPhi");
     Float_t* jetPt    = data.GetPtrFloat("AK5jetPt");
     Float_t* jetEta   = data.GetPtrFloat("AK5jetEta");
     Float_t* jetPhi   = data.GetPtrFloat("AK5jetPhi");
     Float_t* jetMass   = data.GetPtrFloat("AK5jetMass");
   
     vector<bool> over;
     for (int i = 0; i <nJet; i++) {
         over.push_back(1);
     }

     bool checkId =0;
     for (int i = 0; i <nGen; i++) {
         if (genId[i]==15||genId[i]==-15)checkId=1;
     }
     if (checkId)continue;


     for (int i = 0; i <nEl; i++) {
         TLorentzVector ele1;
         ele1.SetPtEtaPhiM(elePt[i], eleEta[i], elePhi[i], 0.000511);
         if(elePt[i]<20)continue;
         for (int j = i + 1; j < nEl; j++) {
             if(elePt[j]<20)continue;
	     TLorentzVector ele2;
	     ele2.SetPtEtaPhiM(elePt[j], eleEta[j], elePhi[j], 0.000511);
	     TLorentzVector  Z = ele1 + ele2; 
	     if ((71<Z.M())&&(Z.M()<111)){
                flag=1;
                q=Z;
	        //q.SetPtEtaPhiM(elePt[i]+elePt[j],eleEta[i]+ eleEta[j], elePhi[i]+elePhi[j], 0.000511*2);
                //cout<<"qset="<<q.Px() <<","<<q.Py() <<endl;
                //cout<<"qvector"<<elePt[i]*cos(elePhi[i])+elePt[j]*cos(elePhi[j])<<","<<elePt[i]*sin(elePhi[i])+elePt[j]*sin(elePhi[j]) <<endl;
                for (int m = 0; m <nJet; m++) {
		    if (((eleEta[i]-jetEta[m])* (eleEta[i]-jetEta[m])+(elePhi[i]-jetPhi[m])*(elePhi[i]-jetPhi[m]))<0.16)over[m]=0;                          if (((eleEta[j]-jetEta[m])* (eleEta[j]-jetEta[m])+(elePhi[j]-jetPhi[m])*(elePhi[j]-jetPhi[m]))<0.16)over[m]=0;
	        }
		break;
             }
          
         }
         if(flag)break;
     }

     if (!flag){ 
         for (int i = 0; i <nMu; i++) {
             TLorentzVector mu1 ; 
	     mu1.SetPtEtaPhiM(muPt[i], muEta[i], muPhi[i], 0.1057);
  	     if(muPt[i]<20)continue;
             for (int j = i + 1; j < nMu; j++) {
                  if(muPt[j]<20)continue;
      	          TLorentzVector mu2;
                  mu2.SetPtEtaPhiM(muPt[j], muEta[j], muPhi[j], 0.1057);
	          TLorentzVector  Z = mu1 + mu2;        
    	          if ((71<Z.M())&&(Z.M()<111)){	  
                      flag=1;
                      q=Z;
                      //q.SetPtEtaPhiM(muPt[i]+muPt[j], muEta[i]+muEta[j],muPhi[i]+ muPhi[j], 0.1057*2);
                      for (int m = 0; m <nJet; m++) {
			if (((muEta[i]-jetEta[m])* (muEta[i]-jetEta[m])+(muPhi[i]-jetPhi[m])*(muPhi[i]-jetPhi[m]))<0.16)over[m]=0;
			if (((muEta[j]-jetEta[m])* (muEta[j]-jetEta[m])+(muPhi[j]-jetPhi[m])*(muPhi[j]-jetPhi[m]))<0.16)over[m]=0;
                      } 
                      break;
                  }
             }
             if(flag)break;
         }
     }

     if (flag){
        
        for (int i = 0; i <nJet; i++) {
           TLorentzVector jet;
           jet.SetPtEtaPhiM(jetPt[i], jetEta[i], jetPhi[i], jetMass[i]);
           if (over[i])u=u+jet;
        }

        qT->Fill(q.Pt());
        u=u+q;
        all->Fill(u.Pt());
        u=u-q;
        uT->Fill(u.Pt());
        double temp=q.Px()*u.Px()+q.Py()*u.Py();
        temp=temp/q.Pt();
        double tempUx=u.Px()-(temp*q.Px())/q.Pt();
        double tempUy=u.Py()-(temp*q.Py())/q.Pt();
        //cout<<"q=("<<q.Px() <<"," << q.Py()<<")"<<endl;
        //cout<<"u=("<<u.Px() <<"," << u.Py()<<")"<<endl;
        //cout<<"up=("<<(temp*q.Px())/q.Pt() <<"," << (temp*q.Py())/q.Pt()<<")"<<endl;
        temp=temp+q.Pt();
        uP->Fill(temp);
        //cout<<"up"<<temp<<endl;
        //temp=u.Pt()-(((q.Px()*u.Px()+q.Py()*u.Py())/q.Pt());
        temp=sqrt(tempUx*tempUx+tempUy*tempUy);
        temp= tempUx>0? temp:  -temp;
        
        //cout<<"upp="<<temp<<endl;
        uPp->Fill(temp);


     }
  }

  TFile* outFile = new TFile("met.root","recreate");    
  qT->Write();
  uT->Write();
  all->Write();
  uP->Write();
  uPp->Write();
  outFile->Close();

}
