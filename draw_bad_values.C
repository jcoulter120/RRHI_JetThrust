// Jennifer Coulter
// July 31st 2015
// Rutgers University, jennifer.coulter@cern.ch
//
// Macro to learn more about bad values of thrust in order to debug thrust_HiForest.C appropriately. 
//

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
#include "TVector3.h"
#include <TROOT.h>

using namespace std;

void draw_bad_values(){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  bool debug = true; 

  //define trees and file
  TFile * fin = TFile::Open("bad_file_100.root");
  TH1F * h_TBadpT = (TH1F*)fin->Get("TBadpT"); 
  TH1F * h_TmajBadpT = (TH1F*)fin->Get("TmajBadpT");
  TH1F * h_TminBadpT = (TH1F*)fin->Get("TminBadpT");
  TH1F * h_TBad = (TH1F*)fin->Get("thrust_bad");
  TH1F * h_TmajBad = (TH1F*)fin->Get("thrust_maj_bad");
  TH1F * h_TminBad = (TH1F*)fin->Get("thrust_min_bad");
  TH1F * h_eta = (TH1F*)fin->Get("etaBad");
  TH1F * h_phi = (TH1F*)fin->Get("phiBad");
  TH1F * h_weight = (TH1F*)fin->Get("weightingBad");
  TH1F * h_pthat = (TH1F*)fin->Get("pthatBad");
  TH1F * h_nref = (TH1F*)fin->Get("nrefBad");
  TH1F * h_jetCount= (TH1F*)fin->Get("jetCountBad");

  TCanvas * canvas = new TCanvas("c","Thrust Test", 1200, 1200);
  //TCanvas * canvas2 = new TCanvas("c2", "Thrust Test", 1200, 1200);
  canvas->Divide(2,3); 

  if (debug) {

    cout << "Thrust Entries: " << endl;
    h_TBad->Print("base");

    cout << "pT Entries: " << endl;
    h_TBadpT->Print("base"); 
  }
 
   //plot with scaling, error bars, and marker styles
  canvas->cd(1)->SetLogy();
  h_TBad->SetTitle("Preliminary Thrust vs. Counts"); 
  h_TBad->SetXTitle("Thrust");
  h_TBad->SetYTitle("dN/dT (1/N)");
  //h_T->SetAxisRange(1e-5,1,"Y");
  h_TBad->GetXaxis()->CenterTitle();
  h_TBad->GetYaxis()->CenterTitle();
  //h_TBad->SetMarkerStyle(20);   h_TmajBad->SetMarkerStyle(21);   h_TminBad->SetMarkerStyle(22);
  h_TBad->SetMarkerStyle(2);   h_TmajBad->SetMarkerStyle(5);   h_TminBad->SetMarkerStyle(3);
  TLegend * a = new TLegend(0.2,.70,.4,.85);
  a->AddEntry(h_TBad, "Thrust", "p");
  a->AddEntry(h_TmajBad, "Thrust Major", "p");
  a->AddEntry(h_TminBad, "Thrust Minor", "p");
  h_TBad->Draw("p");
  h_TminBad->Draw("p&same");
  h_TmajBad->Draw("p&same");
  a->Draw("same");

     //plot with scaling, error bars, and marker styles
  canvas->cd(2)->SetLogy();
  h_TBadpT->SetTitle("pT vs. Counts"); 
  h_TBadpT->SetXTitle("pT");
  h_TBadpT->SetYTitle("Counts");
  //h_T->SetAxisRange(1e-5,1,"Y");
  h_TBadpT->GetXaxis()->CenterTitle();
  h_TBadpT->GetYaxis()->CenterTitle();
  //h_TBad->SetMarkerStyle(20);   h_TmajBad->SetMarkerStyle(21);   h_TminBad->SetMarkerStyle(22);
  h_TBadpT->SetMarkerStyle(2);   h_TmajBadpT->SetMarkerStyle(5);   h_TminBadpT->SetMarkerStyle(3);
  TLegend * v = new TLegend(0.2,.70,.4,.85);
  v->AddEntry(h_TBadpT, "Thrust", "p");
  v->AddEntry(h_TmajBadpT, "Thrust Major", "p");
  v->AddEntry(h_TminBadpT, "Thrust Minor", "p");
  h_TBadpT->Draw("p");
  h_TminBadpT->Draw("p&same");
  h_TmajBadpT->Draw("p&same");
  v->Draw("same");

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
  h_nref->SetAxisRange(1e-5,1000000000,"Y");
  h_nref->Draw();
  h_jetCount->Draw("SAME");			     
  c->Draw("SAME");

  canvas->cd(6)->SetLogy();
  h_pthat->SetTitle("pThat"); 
  h_pthat->SetXTitle("pThat");
  h_pthat->SetYTitle("Counts");
  h_pthat->GetXaxis()->CenterTitle();
  h_pthat->GetYaxis()->CenterTitle();
  h_pthat->Draw("p");

  canvas->Print("bad.pdf");

}
