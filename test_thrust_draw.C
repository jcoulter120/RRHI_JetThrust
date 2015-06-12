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
  
  //define trees and file
  TFile * fin = TFile::Open("test_pp_thrust.root");
  TH1F * h_T = (TH1F*)fin->Get("thrust_scaled");
  TH1F * h_Tmaj = (TH1F*)fin->Get("thrust_maj_scaled");
  TH1F * h_Tmin = (TH1F*)fin->Get("thrust_min_scaled");
  TH1F * h_40 = (TH1F*)fin->Get("thrust_40");
  TH1F * h_60 = (TH1F*)fin->Get("thrust_60");
  TH1F * h_80 = (TH1F*)fin->Get("thrust_80");
  TH1F * h_nref = (TH1F*)fin->Get("nref");
  TH1F * h_jetCount = (TH1F*)fin->Get("jetCount");
  TH1F * h_eta = (TH1F*)fin->Get("eta");
  TH1F * h_phi = (TH1F*)fin->Get("phi");

  TCanvas * c = new TCanvas("c","Thrust Test", 1200, 1200);
  c->Divide(2,3);

  //plot with scaling and marker styles
  c->cd(1)->SetLogy();
  h_T->SetTitle("Preliminary Thrust vs. log(Count)"); 
  h_T->SetXTitle("Thrust");
  h_T->SetYTitle("Counts");
  h_T->GetXaxis()->CenterTitle();
  h_T->GetYaxis()->CenterTitle();
  h_T->SetMarkerStyle(2);   h_Tmaj->SetMarkerStyle(5);   h_Tmin->SetMarkerStyle(3);
  TLegend*legend = new TLegend(0.2,.70,.4,.85);
  legend->AddEntry(h_T, "Thrust", "p");
  legend->AddEntry(h_Tmaj, "Thrust Major", "p");
  legend->AddEntry(h_Tmin, "Thrust Minor", "p");
  h_T->Draw("p");
  h_Tmin->Draw("p&same");
  h_Tmaj->Draw("p&same");
  legend->Draw("same");
  
  //eta bias check plot
  c->cd(2);
  h_eta->SetTitle("Eta vs. dN/dT"); 
  h_eta->SetXTitle("Eta");
  h_eta->SetYTitle("Counts)");
  h_eta->GetXaxis()->CenterTitle();
  h_eta->GetYaxis()->CenterTitle();
  h_eta->Draw(); 

  //phi bias check plot
  c->cd(3);
  h_phi->SetTitle("Phi vs. dN/dT"); 
  h_phi->SetXTitle("Phi (radians)");
  h_phi->SetYTitle("Counts");
  h_phi->GetXaxis()->CenterTitle();
  h_phi->GetYaxis()->CenterTitle();
  h_phi->Draw();
  
  //counts of jets
  c->cd(4)->SetLogy();
  TLegend*g = new TLegend(0.55,.75,.85,.85);
  h_jetCount->SetLineColor(2);   g->AddEntry(h_jetCount,"selected jet count","l");
  h_nref->SetLineColor(4);       g->AddEntry(h_nref,"nref","l");
  h_nref->SetTitle("Preliminary Thrust vs. Log(Counts) for jtpt Cuts"); 
  h_nref->SetXTitle("Number of Jets");
  h_nref->SetYTitle("log(Counts)");
  h_nref->GetXaxis()->CenterTitle();
  h_nref->GetYaxis()->CenterTitle();
  h_nref->Draw();
  h_jetCount->Draw("SAME");			     
  g->Draw("SAME");
  
  // Plots of thrust given Tmaj, Tmin, and T
  c->cd(5)->SetLogy();
  TLegend*p = new TLegend(0.2,.7,.3,.85);
  h_80->SetLineColor(2);   p->AddEntry(h_80,"jtpt80","l");
  h_60->SetLineColor(3);   p->AddEntry(h_60,"jtpt60","l");
  h_40->SetLineColor(4);   p->AddEntry(h_40,"jtpt40","l");
  h_80->SetTitle("Preliminary Thrust vs. Log(Counts) for jtpt Cuts 40, 60, 80"); 
  h_80->SetXTitle("Thrust");
  h_80->SetYTitle("log(Counts)");
  h_80->GetXaxis()->CenterTitle();
  //h_80->SetAxisRange(0, 1000000, "Y");
  h_80->GetYaxis()->CenterTitle();
  h_80->Draw();
  p->Draw("SAME"); 
  h_60->Draw("SAME"); 
  h_40->Draw("SAME");
  c->SaveAs("test_pp_thrust.pdf","RECREATE"); 

}//end of plot thrust




