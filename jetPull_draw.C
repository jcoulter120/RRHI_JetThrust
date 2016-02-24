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

//plot jetPull
void jetPull_draw(){

  gStyle->SetOptStat(0);

  bool debug = true; 
  Int_t jobNum = 60;
  
  //define trees and file
  TFile * fin = TFile::Open(Form("pythia_jetPull_%d.root",jobNum));

  TH1F * h_jetPull = (TH1F*)fin->Get("jetPull");
  TH1F * h_jetPull_2d = (TH1F*)fin->Get("jetPull_2d");
  TCanvas * canvas = new TCanvas("c","jetPull Test", 1200, 1200);
  canvas->Divide(2,3);

  if (debug) {

    cout << "JetPull Entries: " << endl;
    h_jetPull->Print("base");
  }
  
  //plot with scaling, error bars, and marker styles
  canvas->cd(1)->SetLogy();
  h_jetPull->SetTitle("JetPull Magnitude vs. Counts"); 
  h_jetPull->SetXTitle("JetPull");
  h_jetPull->SetYTitle("Counts");
  h_jetPull->SetAxisRange(-5.5,5.5,"X");
  h_jetPull->GetXaxis()->CenterTitle();
  h_jetPull->GetYaxis()->CenterTitle();
  h_jetPull->Draw();

  canvas->cd(2);
  h_jetPull_2d->SetTitle("JetPull"); 
  h_jetPull_2d->SetXTitle("y");
  h_jetPull_2d->SetYTitle("phi");
  h_jetPull_2d->SetAxisRange(-5.5,5.5,"X");
  h_jetPull_2d->GetXaxis()->CenterTitle();
  h_jetPull_2d->GetYaxis()->CenterTitle();
  h_jetPull_2d->Draw("colz");

  canvas->Print(Form("test_pp_jetPull_%d.pdf", jobNum)); 

}//end of plot thrust
