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
#include <TMarker.h>
#include "passMuonID.C"
#include "passElectronID.C"




void met(std::string fin){
  std::vector<string> infiles;
  TString outputFileName;
  bool readOneFile=true;
  if(fin.find(".txt")!= std::string::npos)
    {
      readOneFile=false;
      TString endfix=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"${test%%.txt*}\"",fin.data()));
      outputFileName = Form("met_%s.root",endfix.Data());
      cout << "Output file = " << outputFileName << endl;
    }
  
  if(readOneFile){
      infiles.push_back(fin);
      TString endfix=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"${test%%.root*}\"",fin.data()));
      outputFileName = Form("met_%s.root",endfix.Data());
      cout << "Output file = " << outputFileName << endl;
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

  TH1D* qT = new TH1D("qT", "Removed Vector Boson", 50, 0, 200);
  TH1D* uT = new TH1D("uT", "Sum of All Particels except Removed Vector Boson ", 50, 0, 200);
  TH1D* all = new TH1D("all", "Sum of All Particels ", 50, 0, 200);
  TH1D* uP = new  TH1D("uParaAddqT", "Additon of uT parallel and qT", 50, -200, 200);
  TH1D* uPp = new  TH1D("uPerp", "uT perpendicular ", 50, -200, 200);


  TreeReader data(infiles); // v5.3.12
  data.Print();
  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {


     // for (Long64_t ev = 0; ev < 100; ev++) {
     // print progress
     if (ev % 90000 == 0){
        fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());
     }
     data.GetEntry(ev);

     //----------------------------------------------------------
     
     string muTrig=Form("%s",outputFileName.Data());
     if ( muTrig.find("DoubleMu") != std::string::npos ){

      std::string* trigName = data.GetPtrString("hlt_trigName");
      Int_t* trigResult = data.GetPtrInt("hlt_trigResult");
      const Int_t nsize = data.GetPtrStringSize();

      Bool_t passTrigger = false;

      for(int it = 0; it < nsize; it++){

          std::string thisTrig = trigName[it];
          Int_t results = trigResult[it];

    // muon channel
          if( thisTrig.find("HLT_Mu22_TkMu8") != std::string::npos && results == 1 ){

             passTrigger = true;
             break;

          }

       }
   
       if( !passTrigger ) continue;

  }
     
     //------------------------------------------------



   
     TLorentzVector u;
     u.SetPtEtaPhiM(0,0,0,0);
     TLorentzVector q;
     q.SetPtEtaPhiM(0,0,0,0);

     bool flag=0;

     Int_t nJet= data.GetInt("AK5nJet"); 
     Int_t nMu= data.GetInt("nMu"); 
     Int_t nEl= data.GetInt("nEle");
     
     
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
     //--------------check there is no tou
     if ( muTrig.find("data") ==  std::string::npos ){
         Int_t nGen =data.GetInt("nGenPar");
         Int_t*  genId =data.GetPtrInt("genParId");
         bool checkId =0;
         for (int i = 0; i <nGen; i++) {
             if (genId[i]==15||genId[i]==-15)checkId=1;
         }
         if (checkId)continue;
     }
     //---------------------------------
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
	        //flag=1; //consider electron
                // q=Z;   //consider electron
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
         /* 
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
         */
       Int_t stMu,ndMu;
       if (passMuonID(data,&stMu,&ndMu)){
	 TLorentzVector mu1,mu2,Z; 
	   mu1.SetPtEtaPhiM(muPt[stMu], muEta[stMu], muPhi[stMu], 0.1057);
           mu2.SetPtEtaPhiM(muPt[ndMu], muEta[ndMu], muPhi[ndMu], 0.1057);
           Z=mu1+mu2;
           
	   if(Z.Pt()<80)continue;
          
           if ((71<Z.M())&&(Z.M()<111)){	  
                   flag=1;
                   q=Z;
                   for (int m = 0; m <nJet; m++) {
                       if (((muEta[stMu]-jetEta[m])* (muEta[stMu]-jetEta[m])+(muPhi[stMu]-jetPhi[m])*(muPhi[stMu]-jetPhi[m]))<0.16)over[m]=0;
		       if (((muEta[ndMu]-jetEta[m])* (muEta[ndMu]-jetEta[m])+(muPhi[ndMu]-jetPhi[m])*(muPhi[ndMu]-jetPhi[m]))<0.16)over[m]=0;
                   }
           }
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

  TFile* outFile = new TFile(outputFileName,"recreate");    
  qT->Write();
  uT->Write();
  all->Write();
  uP->Write();
  uPp->Write();
  outFile->Close();

}
