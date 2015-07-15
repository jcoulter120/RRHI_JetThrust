// Jennifer Coulter
// June 9th 2015
// Rutgers University, jennifer.coulter@cern.ch
//
// Plotting macro for sample thrust variables. 
//
//

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
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
#include "TVector3.h"

using namespace std;

//plot thrust
void test_thrust_draw(){
  
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  bool debug = true; 
  
  //define trees and file
  TFile * fin = TFile::Open("test_pp_thrust.root");
  TH1F * h_pT = (TH1F*)fin->Get("pTcut"); 
  TH1F * h_T = (TH1F*)fin->Get("thrust_scaled");
  TH1F * h_Tmaj = (TH1F*)fin->Get("thrust_maj_scaled");
  TH1F * h_Tmin = (TH1F*)fin->Get("thrust_min_scaled");
  TH1F * h_40 = (TH1F*)fin->Get("thrust_40_new");
  TH1F * h_60 = (TH1F*)fin->Get("thrust_60_new");
  TH1F * h_80 = (TH1F*)fin->Get("thrust_80_new");
  TH1F * h_nref = (TH1F*)fin->Get("nref");
  TH1F * h_jetCount = (TH1F*)fin->Get("jetCount");
  TH1F * h_eta = (TH1F*)fin->Get("eta");
  TH1F * h_phi = (TH1F*)fin->Get("phi");
  TH1F * h_pthat = (TH1F*)fin->Get("weighting"); 

  TCanvas * canvas = new TCanvas("c","Thrust Test", 1200, 1200);
  TCanvas * canvas2 = new TCanvas("c2", "Thrust Test", 1200, 1200);
  canvas->Divide(2,3);
  canvas2->Divide(2,3); 

  if (debug) {

    cout << "Thrust Entries: " << endl;
    h_T->Print("base");
    
    cout << "Thrust Maj Entries: " << endl; 
    h_Tmaj->Print("base");
    
    cout << "Thrust Min Entries: " << endl; 
    h_Tmin->Print("base");
    
    cout << "Phi Entries: " << endl; 
    h_phi->Print("base");
    
    cout << "Eta Entries: " << endl;
    h_eta->Print("base");

    cout << "pThat Entries: " << endl;
    h_pthat->Print("base"); 
  }
  
  //plot with scaling, error bars, and marker styles
  canvas->cd(1)->SetLogy();
  h_T->SetTitle("Preliminary Thrust vs. Counts"); 
  h_T->SetXTitle("Thrust");
  h_T->SetYTitle("dN/dT (1/N)");
  //h_T->SetAxisRange(1e-5,1,"Y");
  h_T->GetXaxis()->CenterTitle();
  h_T->GetYaxis()->CenterTitle();
  h_T->SetMarkerStyle(20);   h_Tmaj->SetMarkerStyle(21);   h_Tmin->SetMarkerStyle(22);
  TLegend * a = new TLegend(0.2,.70,.4,.85);
  a->AddEntry(h_T, "Thrust", "p");
  a->AddEntry(h_Tmaj, "Thrust Major", "p");
  a->AddEntry(h_Tmin, "Thrust Minor", "p");
  h_T->Draw("p");
  h_Tmin->Draw("p&same");
  h_Tmaj->Draw("p&same");
  a->Draw("same");
  
  //plot with scaling, error bars, and matching marker styles
  canvas->cd(2)->SetLogy();
  h_pT->SetTitle("Preliminary Thrust vs. Counts"); 
  h_pT->SetXTitle("Thrust");
  h_pT->SetYTitle("(1/N)");
  //h_T->SetAxisRange(1e-5,1,"Y");
  h_pT->GetXaxis()->CenterTitle();
  h_pT->GetYaxis()->CenterTitle();
  //h_T->SetMarkerStyle(2);   h_Tmaj->SetMarkerStyle(5);   h_Tmin->SetMarkerStyle(3);
  //TLegend * b = new TLegend(0.2,.70,.4,.85);
  //b->AddEntry(h_notScaled, "Thrust", "p");
  //b->AddEntry(h_Tmaj, "Thrust Major", "p");
  //b->AddEntry(h_Tmin, "Thrust Minor", "p");
  h_pT->Draw("p");
  //h_Tmin->Draw("p&same");
  //h_Tmaj->Draw("p&same");
  //b->Draw("same");
  
  //eta bias check plot
  canvas->cd(3);
  h_eta->SetTitle("Eta vs. dN/dT"); 
  h_eta->SetXTitle("Eta");
  h_eta->SetYTitle("Counts");
  //h_eta->SetAxisRange(10,110,"Y");
  h_eta->SetAxisRange(-2.5,2.5,"X");
  h_eta->GetXaxis()->CenterTitle();
  h_eta->GetYaxis()->CenterTitle();
  h_eta->Draw(); 

  //phi bias check plot
  canvas->cd(4);
  h_phi->SetTitle("Phi vs. dN/dT"); 
  h_phi->SetXTitle("Phi (radians)");
  h_phi->SetYTitle("Counts");
  //h_phi->SetAxisRange(10,110,"Y");
  h_phi->SetAxisRange(-3.2,3.2,"X");
  h_phi->GetXaxis()->CenterTitle();
  h_phi->GetYaxis()->CenterTitle();
  h_phi->Draw();
  
  //counts of jets
  canvas->cd(5)->SetLogy();
  TLegend*c = new TLegend(0.55,.75,.85,.85);
  h_jetCount->SetLineColor(2);   c->AddEntry(h_jetCount,"selected jet count","l");
  h_nref->SetLineColor(4);       c->AddEntry(h_nref,"nref","l");
  h_nref->SetTitle("Preliminary Thrust vs. Counts for pT Cuts"); 
  h_nref->SetXTitle("Number of Jets");
  h_nref->SetYTitle("Counts");
  h_nref->GetXaxis()->CenterTitle();
  h_nref->GetYaxis()->CenterTitle();
  h_nref->SetAxisRange(1e-5,10000,"Y");
  h_nref->Draw();
  h_jetCount->Draw("SAME");			     
  c->Draw("SAME");
  
  // Plots of thrust given Tmaj, Tmin, and T
  canvas->cd(6)->SetLogy();
  TLegend*d = new TLegend(0.2,.7,.3,.85);

  h_80->SetLineColor(2);   d->AddEntry(h_80,"jtpt80","l");
  h_60->SetLineColor(3);   d->AddEntry(h_60,"jtpt60","l");
  h_40->SetLineColor(4);   d->AddEntry(h_40,"jtpt40","l");
  h_80->SetTitle("Preliminary Thrust vs. Counts for jtpt Cuts 40, 60, 80"); 
  h_80->SetXTitle("Thrust");
  h_80->SetYTitle("dN/dT (1/N)");
  h_80->GetXaxis()->CenterTitle();
  h_80->GetYaxis()->CenterTitle();
  h_80->SetAxisRange(1e-5,1,"Y");
  h_80->Draw("p");
  d->Draw("SAME"); 
  h_60->Draw("p&SAME"); 
  h_40->Draw("p&SAME");

  canvas2->cd(1)->SetLogy();
  h_pT->SetTitle("pT after Cuts"); 
  h_pT->SetXTitle("pT");
  h_pT->SetYTitle("Counts");
  h_pT->GetXaxis()->CenterTitle();
  h_pT->GetYaxis()->CenterTitle();
  h_pT->Draw("p");

  canvas2->cd(2)->SetLogy();
  h_pthat->SetTitle("pThat"); 
  h_pthat->SetXTitle("pThat");
  h_pthat->SetYTitle("Counts");
  h_pthat->GetXaxis()->CenterTitle();
  h_pthat->GetYaxis()->CenterTitle();
  h_pthat->Draw("p");

  canvas2->Print("pTcut.png"); 
  canvas->Print("test_pp_thrust.pdf"); 

}//end of plot thrust




