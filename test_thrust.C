// Jennifer Coulter
// May 22th 2015
// Rutgers, jennifer.coulter@cern.ch
//
// Test macro for plotting thrust, an event shape variable.
//
// May 25th -> added for loops for events and jets
// May 26th -> added plot with respect to dN/dT 

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


using namespace std;
void test_thrust(){

  bool debug = false;
  
  //define trees and file
  TFile * fin = TFile::Open("pp_2013_data_testfile.root");
  TFile * save_File = new TFile("test_pp_thrust.root","RECREATE");
  TTree * t = (TTree*)fin->Get("ak3PFJetAnalyzer/t");
  TCanvas * c = new TCanvas("c","Thrust Test", 1200, 1200);
  c->Divide(2,1);
  TH1F * h_thrust = new TH1F("thrust", "", 50,0,1);

  Double_t px[1000];
  Double_t py[1000];
  Double_t pz[1000];
  Float_t pt[1000];
  Float_t eta[1000];
  Float_t phi[1000];
  Float_t theta[1000];
  Float_t thrust[50000]; 
  Int_t nref;
  Int_t p_tot = 0;
  Float_t dot = 0;
  Double_t mag = 0;
  Double_t axis_theta = 0;
  Double_t thrust_temp = 0;
  Double_t thrust_max = 0;
  Double_t nT_mag = 0; 
  
  //Set branches of the tree 
  t->SetBranchAddress("jtpt", &pt);
  t->SetBranchAddress("jteta", &eta);
  t->SetBranchAddress("jtphi", &phi);
  t->SetBranchAddress("nref", &nref);
  
  Long64_t nentries = t->GetEntries();
  if(debug) nentries = 1000;
  
  //event loop
for(Long64_t nentry = 0; nentry<nentries; ++nentry){
    
    t->GetEvent(nentry);

    //make a selection cut???
    bool select = false;
    for(int k = 0; k < nref; k++){

      if(pt[k] > 30){ select = true;}
      
    }

    if(!select) continue;
    
    //max axis loop
    for(Long64_t naxis = 0; naxis < nref; ++naxis){
      
      thrust_max = 0; 
      
      //calculates theta for this jet
      axis_theta = 2*TMath::ATan(exp(-1*eta[naxis]));
      if(debug) cout<<"axis theta  = " << axis_theta << endl; 
      
      //calculates axis for this particular jet
      TVector3 nT (TMath::Sin(axis_theta) * TMath::Cos(phi[naxis]), TMath::Sin(phi[naxis]) * TMath::Sin(axis_theta), TMath::Cos(phi[naxis]));
      
      //normalize the thrust axis
      nT_mag = TMath::Sqrt(nT(0)*nT(0) + nT(1)*nT(1) + nT(2)*nT(2)); 
      nT(0) = nT(0)/nT_mag;
      nT(1) = nT(1)/nT_mag;
      nT(2) = nT(2)/nT_mag;
      
      if(debug) cout<<"axis  = " << nT(0) << " " << nT(1) << " " << nT(2) << endl;
      if(debug) cout<<"axis mag  = " << nT.Mag() << endl;

      //resets for next jet loop
      p_tot = 0;
      dot = 0;
      mag = 0; 
      
      //jet loop
      for(Long64_t njet = 0; njet < nref; ++njet){
	
	//calculate theta from eta
	theta[njet] = 2*TMath::ATan(exp(-1*eta[njet]));
	//cout<<thrust_max<<endl;	
	thrust_temp = 0;
	
	//calculate px, py, pz
	px[njet] = pt[njet]*TMath::Cos(phi[njet]);
	if(debug) cout<<"px = " << px[njet] << endl; 
	py[njet] = pt[njet]*TMath::Sin(phi[njet]);
	if(debug) cout<<"py = " << py[njet] << endl; 
	pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	if(debug) cout<<"pz = " << pz[njet] << endl;
	
	//define momentum three vector
	TVector3 p3 (px[njet], py[njet], pz[njet]);
	
	//sum the total p from the individual p magnitudes
	mag += p3.Mag();
	if(debug) cout<<"mag = " << mag << endl;
	
	//dots the two vectors
	dot += TMath::Abs(px[njet]*nT(0) + py[njet]*nT(1) + pz[njet]*nT(2));
	//dot += (px[njet]*nT(0) + py[njet]*nT(1) + pz[njet]*nT(2));	
	if(debug) cout<<"dot = " << dot << endl;	
	
      }//end jet loop
      
      //calculate the thrust
      thrust_temp = ((dot)/mag);

      //Compare to see if this axis is a new maximum 
      if(debug) cout<< "temp thrust = " << thrust_temp << endl; 
      if(debug) cout<< "max thrust = " << thrust_max << endl;
      
      if(thrust_temp>thrust_max){
	thrust_max = thrust_temp;
	if(debug) cout<< " thrust max "<< thrust_max << endl; 
      }
      
    }//end axis loop
    
    h_thrust->Fill(thrust_max);
    if(debug) cout<< "entry value = " << thrust_max;  
    
  }//end of event loop
  
  //Set up histogram
  c->cd(1)->SetLogy();
  h_thrust->SetTitle("Preliminary Thrust vs. log(Count)"); 
  h_thrust->SetXTitle("Thrust");
  h_thrust->SetYTitle("log(Counts)");

  h_thrust->Print("base");
  cout<<"histogram mean = "<<h_thrust->GetMean()<<endl;
  h_thrust->Draw();

  //Create the plot for Thrust vs. dN/dT
  //define a second histogram
  TH1F * h_T = new TH1F("thrust1", "", 50,0,1);
  Float_t bin = h_T->GetBinWidth(0); 
  
  //h_thrust->Sumw2();
  for (int i=0;i<h_thrust->GetNbinsX();i++)
    {
      Float_t val = (h_thrust->GetBinContent(i+1)) - (h_thrust->GetBinContent(i));
      
      //Float_t valErr = h_thrust->GetBinError(i);
      val = val/bin; 
      //valErr/=h_thrust->GetBinWidth(i);
      // h_T->SetBinContent(i,val);
      for(int p = 0; p < val; p++){

      	h_T->AddBinContent(i); 
      }
      //cout<<"val = " << val << endl;
      //h_T->SetBinError(i,valErr);
    }

  h_T->GetXaxis()->CenterTitle();
  h_T->GetYaxis()->CenterTitle();

  //Set up second histogram
  c->cd(2)->SetLogy();
  h_T->SetTitle("Preliminary Thrust vs. dN/dT"); 
  h_T->SetXTitle("Thrust");
  h_T->SetYTitle("dN/dT");

  h_T->Print("base");
  cout<<"histogram mean = "<<h_T->GetMean()<<endl;
  h_T->Draw();

  save_File->Write();
  c->SaveAs("test_pp_thrust.root","RECREATE");
  
}//end of macro
