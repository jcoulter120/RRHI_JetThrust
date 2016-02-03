// Jennifer Coulter
// May 19th 2015
// Rutgers, jbc120@scarletmail.rutgers.edu
//
//
// test macro to read hiForest and plot jet variables
//
//May 19th 2015 --> added for loop rendition of histogram plotting

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

using namespace std;

int whichBin(Float_t pthat){
  
  int bin; 
  if(pthat < 30) bin = 0;
  else if(pthat < 50) bin = 1;
  else if(pthat < 80) bin = 2;
  else if(pthat < 120) bin = 3;
  else if(pthat < 170) bin = 4;
  else if(pthat < 220) bin = 5;
  else if(pthat < 280) bin = 6;
  else if(pthat < 330) bin = 7;
  else if(pthat < 400) bin = 8;
  else if(pthat < 460) bin = 9;
  else if(pthat < 540) bin = 10;
  else{ bin = 11; }
  return bin; 
}

void test_macro(){
  
  /* COUNT UP FILE NAMES ================================= */
  //string input_file = "Text_Files/jewel_test_med10_5k.txt";
  string input_file = "Text_Files/jewel_files_1k.txt";
  
  //DEFINE VARIABLES
  TFile * file;
  Int_t startfile = 0;
  Int_t endfile = 100;
  ifstream count(input_file.c_str(), ifstream::in);
  TTree * t;
  TTree * nt;
  Float_t pThat; 
  Int_t fileCount = 0;
  string * filename = new string[1350];
  int * eventTally = new int[1350];
  TH1F * hJet_pT = new TH1F("hJet_pT", "", 100,0,1000);
  TH1F * h_nref = new TH1F("nref", "", 12, 0, 12);
  TH1F * h_weight = new TH1F("weighting", "", 1200, 0, 1200);
  TH1F * h_unweight = new TH1F("unweighted", "", 1200, 0, 1200);
  Float_t pt[1000];
  Int_t nref;
  bool debug = false;
  Double_t weightSTAT;
  double pthatBinning[] = {15,30,50,80,120,170,220,280,330,400,460,540};
  //double xsSTAT[] = {286240.0,444780.0,494830.0,514200.0,554520.0,482280.0,480060.0,309090.0,0,0,114700.0,53175.0};
  //double xsSTAT[] = {4878812.4,4218696.2,4659866.4,4049616.0,2934218.1,1695768.4,1302251.3,718081.6,0,0,168207.79,71406.91};
  double xsSTAT[] = {4878812.4,4218696.2,4659866.4,4049616.0,2934218.1,1695768.4,1302251.3,718081.6,578537.1,291444.26,168207.79,71406.91};
  
  
  //count up the total number of files and save their names to an array of filenames for initialization purposes
  
  string line;
  while(getline(count, line)){
    filename[fileCount] = line;
    if (debug) cout << filename[fileCount] << endl; 
    fileCount++;
  }
  count.close();
  /* ===================================================== */
  
  //OPEN FILE
  for(int ifile = startfile; ifile < endfile; ifile++){
    
    string str = filename[ifile];
    file = TFile::Open(str.c_str());
    if(file == NULL) { cout << "File *** " << filename[ifile] << " *** does not exist." << endl;  continue; }
    
    cout << "\n **** =========================== New File ================================= **** \n ";
    cout << "File Name: " << filename[ifile] << endl;
    cout << "File Number: " << ifile << "/" << startfile << " to " << endfile << endl;
    
    //SET BRANCHES
    t = (TTree*)file->Get("dijet/t3");
    nt = (TTree*)file->Get("dijet/nt");
    
    t->SetBranchAddress("jtpt", &pt);
    t->SetBranchAddress("nref", &nref);
    nt->SetBranchAddress("pthat", &pThat);
    
    t->AddFriend(nt);
    
    Long64_t nentries =  t->GetEntries();
    cout<< "Number of events "<<t->GetEntries()<<endl;

    /* COUNT UP THE NUMBER OF EVENTS IN EACH FILE ========= */
    for(Int_t h = 0; h < nentries; ++h){
      t->GetEvent(h);
      if(nref != 0){
	eventTally[ifile]++;
      }
      else{
	//break; 
      }
    }
    /* ==================================================== */
    cout << "Number of Decent Events: " << eventTally[ifile] << endl;
    
    if(debug) nentries = 10;
    
    /* EVENT LOOP ======================================= */
    for(Long64_t nentry = 0; nentry < nentries; ++nentry){
      
      t->GetEvent(nentry);

      //if(nref == 0) {break;}

      weightSTAT = xsSTAT[whichBin(pThat)]/nentries;
      weightSTAT = xsSTAT[whichBin(pThat)]/eventTally[ifile];
      h_nref->Fill(nref);
      h_weight->Fill(pThat,weightSTAT);
      h_unweight->Fill(pThat);

      /* START THE JET LOOP ============================= */ 
      for(int jentry = 0; jentry < nref; ++jentry){
	
	//hJet_pT->Fill(pt[jentry], weightSTAT);
	hJet_pT->Fill(pt[jentry]);
	//if(printDebug) cout<<"pt = "<< pt[jentry] << endl;
	
      }/* END JET LOOP ================================== */
    }/* END EVENT LOOP ================================== */
  }//END OF FILE LOOP
  /* PLOT THE JTPT ==================================== */
  hJet_pT->Print("base");
  cout<<"histogram mean = "<<hJet_pT->GetMean()<<endl;
  
  h_nref->Print("base");
  cout<<"histogram mean = "<<h_nref->GetMean()<<endl;
  
  TCanvas * c2 = new TCanvas("c2", "Simple Jet Variables 2", 1200, 1200);
  c2->Divide(2,2);
  
  c2->cd(1)->SetLogy();
  hJet_pT->SetTitle("Jet pT");
  hJet_pT->SetXTitle("Jet pT (GeV/c)");
  hJet_pT->SetYTitle("Counts");
  hJet_pT->Draw("p");
  hJet_pT->SetMarkerStyle(2);
  
  c2->cd(2)->SetLogy();
  h_nref->SetTitle("nref");
  h_nref->SetXTitle("nref");
  h_nref->SetYTitle("Counts");
  h_nref->Draw("p");
  h_nref->SetMarkerStyle(2);
  
  c2->cd(3)->SetLogy();
  h_unweight->SetTitle("pT without Weighting"); 
  h_unweight->SetXTitle("pThat");
  h_unweight->SetYTitle("Counts");
  h_unweight->SetMarkerStyle(2);
  h_unweight->Draw("p");
  
  c2->cd(4)->SetLogy();
  h_weight->SetTitle("Weighted");
  h_weight->SetXTitle("nref");
  h_weight->SetYTitle("Counts");
  h_weight->Draw("p");
  h_weight->SetMarkerStyle(2);
  
  /* ================================================== */
  
  //new_File->Write();
  //c2->SaveAs("test_pp_jetvariables_with_fors.pdf","RECREATE");
  
}// end of macro

