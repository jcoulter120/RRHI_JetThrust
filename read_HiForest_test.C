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

void read_HiForest_test(Int_t radius = 3, char * algo = (char*)"PF"){

  bool printDebug = false;

  //Define file and tree
  TFile * fin = TFile::Open("pp_2013_data_testfile.root");
  TFile * new_File = new TFile("test_pp_jetVariable_plots.root","RECREATE");
  TTree * jet = (TTree*)fin->Get(Form("ak%d%sJetAnalyzer/t",radius,algo));
   
  //Part 1: Plotting Manually
  
  //Define variables
  TH1F * hpT = new TH1F("hpT","",100,0,1000);
  TH1F * heta = new TH1F("heta","",60,-3,+3);
  TH1F * hphi = new TH1F("hphi","",60,-3.15,+3.15);
  TH2F * h2d = new TH2F("h2d","pT vs eta for pT>20",30,-2,+2,100,0,500);

  //Create canvas
  TCanvas * c1 = new TCanvas("c1","Simple Jet Variables 1",1200,1200);
  c1->Divide(2,2);

  //Plotting
  c1->cd(1);
  c1->cd(1)->SetLogy();
  jet->Draw("jtpt>>hpT");
  hpT->SetTitle("pp ak3PF jet pT");
  hpT->SetXTitle("Jet p_{T} (GeV/c)");
  hpT->SetYTitle("Counts");
  hpT->Draw();

  c1->cd(2);
  jet->Draw("jteta>>heta");
  heta->SetTitle("pp ak3PF jet eta");
  heta->SetXTitle("Jet Eta");
  heta->SetYTitle("Counts");
  heta->Draw();
  
  c1->cd(3);
  jet->Draw("jtphi>>hphi");
  hphi->SetTitle("pp ak3PF jet phi");
  hphi->SetXTitle("Jet Phi (rad)");
  hphi->SetYTitle("Counts");
  hphi->Draw();
  
  c1->cd(4);
  jet->Draw("jtpt:jteta>>h2d","jtpt>20","goff");
  h2d->SetTitle("pp ak3PF jet jteta v. jtpt");
  h2d->SetXTitle("jet Eta");
  h2d->SetYTitle("jet pT (GeV/c)");
  h2d->Draw("colz");

  c1->SaveAs("test_pp_jetvariables.pdf","RECREATE");
  
  //Part 2: Drawing Using Loops

  //Definitions
  TH1F * hJet_pT = new TH1F("hJet_pT", "", 100,0,1000);
  TH1F * hJet_eta = new TH1F("hJet_eta", "", 60,-3,+3);
  TH1F * hJet_phi = new TH1F("hJet_phi", "", 60,-3.15,+3.15);
  TH2F * hJet_2d = new TH2F("hJet_2d","Jet pT vs Eta for pT>20",30,-3,+3,100,0,500);
   
  Float_t pt[1000];
  Float_t eta[1000];
  Float_t phi[1000];
  Int_t nref;
  
  //Set branches of the tree 
  jet->SetBranchAddress("jtpt", &pt);
  jet->SetBranchAddress("jteta", &eta);
  jet->SetBranchAddress("jtphi", &phi);
  jet->SetBranchAddress("nref", &nref);

  Long64_t nentries =  jet->GetEntries();
  cout<< "Number of events "<<jet->GetEntries()<<endl;

  if(printDebug) nentries = 10;
  
  //Start entry loop
  for(Long64_t nentry = 0; nentry < nentries; ++nentry){
    
    jet->GetEvent(nentry);
    //if(nentry%1000 == 0) cout << nentry << "/" << nentries <<endl;
    
    //start the jet loop
    for(int jentry = 0; jentry < nref; ++jentry){

      hJet_pT->Fill(pt[jentry]);
      //if(printDebug) cout<<"pt = "<< pt[jentry] << endl;

      hJet_eta->Fill(eta[jentry]);
      //if(printDebug) cout<<"eta = "<< eta[jentry] << endl;

      hJet_phi->Fill(phi[jentry]);
      //if(printDebug) cout<<"phi = "<< phi[jentry] << endl;

      hJet_2d->Fill(eta[jentry],pt[jentry]);
      //if(printDebug) cout<<"pt = "<< pt[jentry] << endl;
      //if(printDebug) cout<<"eta  = "<< eta[jentry] << endl;

    }//end jet loop
   
  }//end entry loop

  // check histograms (Prints histogram information)
  hJet_2d->Print("base");
  cout<<"histogram mean = "<<hJet_pT->GetMean()<<endl;
  
  // plotting
  TCanvas * c2 = new TCanvas("c2", "Simple Jet Variables 2", 1200, 1200);
  c2->Divide(2,2);
  
  c2->cd(1)->SetLogy();
  hJet_pT->SetTitle("pp ak3PF Jet pT");
  hJet_pT->SetXTitle("Jet pT (GeV/c)");
  hJet_pT->SetYTitle("Counts");
  hJet_pT->Draw();

  c2->cd(2);
  hJet_eta->SetTitle("pp ak3PF Jet Eta");
  hJet_eta->SetXTitle("Jet Eta");
  hJet_eta->SetYTitle("Counts");
  hJet_eta->Draw();
  
  c2->cd(3);
  hJet_phi->SetTitle("pp ak3PF Jet Phi");
  hJet_phi->SetXTitle("Jet Phi (Rad)");
  hJet_phi->SetYTitle("Counts");
  hJet_phi->Draw();

  c2->cd(4);
  hJet_2d->SetTitle("pp ak3PF Jet Eta vs. pT for pT>20");
  hJet_2d->SetXTitle("Jet Eta");
  hJet_2d->SetYTitle("Jet pT (GeV/c)");
  hJet_2d->Draw("colz");

  new_File->Write();
  c2->SaveAs("test_pp_jetvariables_with_fors.pdf","RECREATE");
  
}// end of macro
