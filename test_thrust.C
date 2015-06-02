// Jennifer Coulter
// May 22th 2015
// Rutgers, jennifer.coulter@cern.ch
//
// Test macro for plotting thrust, an event shape variable.
//
// May 25th -> added for loops for events and jets
// May 26th -> added plot with respect to dN/dT
// May 29th -> corrected error in variable reset, identified issue with axes definition
// May 31st -> added thrust min and max plotting

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

//Can I define a plane class....?
/*
TVector3 projection(TVector3 nT, TVector3 jaxis){

  TVector3 proj;

  //Use TVector3 to find an orthogonal vector
  TVector3 v1 = nT.Orthogonal();

  //Find a vector in the same plane that is still ortho to nT
  TVector3 v2 = nT.Cross(v1); 

  //GS Orthogonalize
  TVector3 u1 = v1;
  TVector3 u2 = v2 - ((v2.Dot(u1))/u1.Dot(u1))*u1; 

  //Find the projection of a jet onto this subspace
  proj = jaxis.Dot(u1)/(u1.Dot(1))*(u1) + (jaxis.Dot(u2)/u2.Dot(u2))*u2; 
  
  //return this projection
  return proj; 
}
*/
//creates histograms in terms of Thrust vs. dN/dT
TH1F* rebin(TH1F * hist, const char * name){

  TH1F* h_return = new TH1F(name, "", hist->GetNbinsX(), 0,1);
  Float_t bin = hist->GetBinWidth(1);
  
  for (int i=0;i<=hist->GetNbinsX();++i){
    Float_t val = hist->GetBinContent(i);
    val = val/bin;
    //Float_t valErr = h_return->GetBinError(i);
      //valErr/=h_return->GetBinWidth(i);
      for(int p = 0; p < val; p++){
	h_return->AddBinContent(i);
	//h_return->SetBinError(i,valErr);
    }
  }
  return h_return;
}

using namespace std;

//plot thrust
void test_thrust(){

  bool debug = false;
  
  //define trees and file
  TFile * fin = TFile::Open("pp_2013_data_testfile.root");
  TFile * save_File = new TFile("test_pp_thrust.root","RECREATE");
  TTree * t = (TTree*)fin->Get("ak3PFJetAnalyzer/t");
  TTree * hiEvt = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  TTree * skim = (TTree*)fin->Get("skimanalysis/HltTree"); 
  TCanvas * c = new TCanvas("c","Thrust Test", 1200, 800);
  
  c->Divide(2,1);
  TH1F * h_thrust = new TH1F("thrust", "", 50,0,1);
  TH1F * h_min = new TH1F("thrust min", "", 50,0,1);
  TH1F * h_maj = new TH1F("thrust maj", "", 50,0,1);

  Double_t px[1000];
  Double_t py[1000];
  Double_t pz[1000];
  Float_t pt[1000];
  Float_t eta[1000];
  Float_t phi[1000];
  
  Int_t nref;
  Float_t dot = 0;
  Double_t mag = 0;
  Double_t thrust_temp = 0;
  Double_t thrust_max = 0;
  Double_t nT_mag = 0;
  Double_t dot_maj = 0;
  Double_t dot_min = 0;
  Double_t min_temp = 0;
  Double_t maj_temp = 0;
  Double_t thrust_maj =0;
  Double_t thrust_min = 0;
  Float_t vz;
  Int_t halo;
  Int_t noise; 
  
  //Set branches of the tree 
  t->SetBranchAddress("jtpt", &pt);
  t->SetBranchAddress("jteta", &eta);
  t->SetBranchAddress("jtphi", &phi);
  t->SetBranchAddress("nref", &nref);
  t->AddFriend(hiEvt);
  t->SetBranchAddress("vz", &vz);
  t->AddFriend(skim);
  t->SetBranchAddress("pPAcollisionEventSelectionPA",&halo);
  t->SetBranchAddress("pHBHENoiseFilter", &noise);
  
  Long64_t nentries = t->GetEntries();
  if(debug) nentries = 1000;
  
  //event loop
  for(Long64_t nentry = 0; nentry<nentries; ++nentry){
    
    t->GetEvent(nentry);

    //make selection cuts
    bool select = true;
    bool pTSelect = false;
    bool etaSelect = false;
    
    //pT cut
    for(int k = 0; k < nref; k++){
      if(pt[k] >  30){cout<< "selected pT"<<endl; pTSelect = true;}
      if((TMath::Abs(eta[k])<2)) {cout << "selected eta"<< endl; etaSelect = true; }
    }

    cout<< "noise = " << noise << " halo = " << halo << endl; 
    if((TMath::Abs(vz) > 15)||(noise=0)||(halo=0)) {select = false;}
    
    if(!select||!pTSelect||!etaSelect) continue;
    if(debug) cout<< " \n ******* New Event ******** " << endl;
    thrust_max = 0; 
    
    //max axis loop
    for(Long64_t naxis = 0; naxis < nref; ++naxis){

      if(debug) cout<< " \n --------- New Test Axis --------- " << endl; 
      
      thrust_temp = 0; 
      
      px[naxis] = pt[naxis]*TMath::Cos(phi[naxis]);
      py[naxis] = pt[naxis]*TMath::Sin(phi[naxis]);
      pz[naxis] = pt[naxis]*TMath::SinH(eta[naxis]);
      
      //calculates axis for this particular jet
      TVector3 nT (px[naxis], py[naxis], pz[naxis]);

      //define the plane perpendicular to nT
      
      if(debug) cout<<"Test Axis UNNORM'D = {" << nT(0) << ", " << nT(1) << ", " << nT(2)<< "}" << endl;	    
      if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[naxis]<<"\n \t eta = "<<eta[naxis]<<"\n \t phi = "<<phi[naxis]<<endl;
      
      //normalize the thrust axis
      nT_mag = TMath::Sqrt(nT(0)*nT(0) + nT(1)*nT(1) + nT(2)*nT(2)); 
      nT(0) = nT(0)/nT_mag;
      nT(1) = nT(1)/nT_mag;
      nT(2) = nT(2)/nT_mag;

      TVector3 maj_axis = nT.Orthogonal();
      TVector3 min_axis = nT.Cross(maj_axis);
      
      if(debug) cout<<"Test Axis NORM'D  = {" << nT(0) << ", " << nT(1) << ", " << nT(2)<< "}" << endl;
      
      //resets for next jet loop
      dot = 0;
      dot_maj = 0;
      dot_min = 0;
      mag = 0; 
      
      //jet loop
      for(Long64_t njet = 0; njet < nref; ++njet){
	
	if(debug) cout<< " \n --------- New Jet --------- " << endl; 
	if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[njet]<<"\n \t eta = "<<eta[njet]<<"\n \t phi = "<<phi[njet]<<endl;

	//calculate px, py, pz
	px[njet] = pt[njet]*TMath::Cos(phi[njet]);
	py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
	pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	
	//define momentum three vector
	TVector3 p3 (px[njet], py[njet], pz[njet]);
	if(debug) cout<<"\nJet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;

	//dots the two vectors for Thrust, Tmin and Tmaj
	dot += TMath::Abs(p3.Dot(nT)); 
	if(debug) cout<<"dot sum = " << dot << endl;

	dot_maj += TMath::Abs(p3.Dot(maj_axis)); 
	if(debug) cout<<"dot sum = " << dot << endl;

	dot_min += TMath::Abs(p3.Dot(min_axis)); 
	if(debug) cout<<"dot sum = " << dot << endl;
	
	//sum the total p from the individual p magnitudes
	mag += TMath::Abs(p3.Mag());
	if(debug) cout<<"mag sum = " << mag << endl;
	
      }//end jet loop
      
      //calculate the thrust
      thrust_temp = ((dot)/mag);
      maj_temp = dot_maj/mag;
      min_temp = dot_min/mag;

      //Compare to see if this axis is a new maximum 
      if(debug) cout<< "\ntemp thrust = " << thrust_temp << endl; 
      
      if(thrust_temp>thrust_max){
	thrust_max = thrust_temp;
	thrust_maj = maj_temp;
	thrust_min = min_temp; 
	if(debug) cout<< "max thrust = " << thrust_max << endl;
      }
      
    }//end axis loop
    h_thrust->Fill(thrust_max);
    h_min->Fill(thrust_min);
    h_maj->Fill(thrust_maj);
    //send the max axis to the other two functions
    
    if(thrust_max < .5) cout<< " thrust entry & "<< thrust_max << endl; 
    
  }//end of event loop
  
  //Set up histograms
  c->cd(1)->SetLogy();
  h_thrust->SetTitle("Preliminary Thrust vs. log(Count)"); 
  h_thrust->SetXTitle("Thrust");
  h_thrust->SetYTitle("log(Counts)");
  h_thrust->GetXaxis()->CenterTitle();
  h_thrust->GetYaxis()->CenterTitle();
  h_thrust->SetLineColor(2);
  h_maj->SetLineColor(3);
  h_min->SetLineColor(4);
  TLegend*legend = new TLegend(0.15,.75,.3,.8);
  legend->AddEntry(h_thrust, "Thrust", "l");
  //legend->AddEntry(h_maj, "Thrust Major", "l");
  //legend->AddEntry(h_min, "Thrust Minor", "l");
  
  h_thrust->Print("base");
  cout<<"histogram mean = "<<h_thrust->GetMean()<<endl;
  gStyle->SetOptStat(0);
  h_thrust->Draw();

  h_maj->Print("base");
  cout<<"histogram mean = "<<h_maj->GetMean()<<endl;
  gStyle->SetOptStat(0);
  //h_maj->Draw("SAME");
  
  h_min->Print("base");
  cout<<"histogram mean = "<<h_min->GetMean()<<endl;
  gStyle->SetOptStat(0);
  //h_min->Draw("SAME");
  legend->Draw("SAME");

  /*
  c->cd(3)->SetLogy();
  h_maj->SetTitle("Preliminary Thrust Major vs. log(Count)"); 
  h_maj->SetXTitle("Thrust");
  h_maj->SetYTitle("log(Counts)");
  h_maj->GetXaxis()->CenterTitle();
  h_maj->GetYaxis()->CenterTitle();
  
  
  c->cd(5)->SetLogy();
  h_min->SetTitle("Preliminary Thrust Minor vs. log(Count)"); 
  h_min->SetXTitle("Thrust");
  h_min->SetYTitle("log(Counts)");
  h_min->GetXaxis()->CenterTitle();
  h_min->GetYaxis()->CenterTitle();
  */
  
  //Create the plot for Thrust vs. dN/dT
  //define histograms
  TH1F * h_T = rebin(h_thrust, "thrust");
  TH1F * h_Tmaj = rebin(h_maj, "thrustmaj");
  TH1F * h_Tmin = rebin(h_min, "thrustmin");
  
  //Set up second round of histograms
  c->cd(2)->SetLogy();
  TLegend*leg = new TLegend(0.15,.75,.3,.8);
  leg->AddEntry(h_T,"Thrust","p");
  //leg->AddEntry(h_Tmaj,"Thrust Major","p");
  //leg->AddEntry(h_Tmin,"Thrust Minor","p");
  h_T->SetTitle("Preliminary Thrust vs. dN/dT"); 
  h_T->SetXTitle("Thrust");
  h_T->SetYTitle("dN/dT");
  h_T->GetXaxis()->CenterTitle();
  h_T->GetYaxis()->CenterTitle();
  h_T->SetMarkerStyle(2);
  h_Tmaj->SetMarkerStyle(5);
  h_Tmin->SetMarkerStyle(3);
  h_T->SetOption("P");
  h_Tmaj->SetOption("P");
  h_Tmin->SetOption("P");
  h_T->Draw();
  leg->Draw("SAME");
  h_T->Print("base");
  cout<<"histogram mean = "<<h_T->GetMean()<<endl;
  //h_Tmaj->Draw("SAME");
  //h_Tmin->Draw("SAME");
  
  /*
  //c->cd(4)->SetLogy();
  h_Tmaj->SetTitle("Preliminary Thrust Major vs. dN/dT"); 
  h_Tmaj->SetXTitle("Thrust");
  h_Tmaj->SetYTitle("dN/dT");
  h_Tmaj->GetXaxis()->CenterTitle();
  h_Tmaj->GetYaxis()->CenterTitle();
  //h_T->SetMarkerStyle(21);
  //h_Tmaj->Draw("SAME");
  
  // c->cd(6)->SetLogy();
  h_Tmin->SetTitle("Preliminary Thrust Minor vs. dN/dT"); 
  h_Tmin->SetXTitle("Thrust");
  h_Tmin->SetYTitle("dN/dT");
  h_Tmin->GetXaxis()->CenterTitle();
  h_Tmin->GetYaxis()->CenterTitle();
  //h_T->SetMarkerStyle(21);
  //h_Tmin->Draw("SAME");
  
  */
  
  //check the histograms
  if(debug){

    h_Tmaj->Print("base");
    cout<<"histogram mean = "<<h_Tmaj->GetMean()<<endl;
    h_Tmin->Print("base");
    cout<<"histogram mean = "<<h_Tmin->GetMean()<<endl;
  }
  
  save_File->Write();
  c->SaveAs("test_pp_thrust.pdf","RECREATE");
  c->SaveAs("test_pp_thrust.root","RECREATE");
  
  
}//end of plot thrust






//Things I might need for updates to this macro: 
      
      //calculates theta for this jet
      // axis_theta = 2*TMath::ATan(exp(-1*eta[naxis]));
      //if(debug) cout<<"axis theta  = " << axis_theta << endl;

      //h_T->SetBinError(i,valErr);
      //Float_t valErr = h_thrust->GetBinError(i);
      //valErr/=h_thrust->GetBinWidth(i);
      // h_T->SetBinContent(i,val);

      //TVector3 nT (TMath::Sin(axis_theta) * TMath::Cos(phi[naxis]), TMath::Sin(phi[naxis]) * TMath::Sin(axis_theta), TMath::Cos(phi[naxis]));


//MONDAY TO-DO

//FIX THE LEGENDS
//PUT THEM IN THE POWERPOINT
//FIX BIN SIZE
//MAKE ALL THE SECOND PLOT OUT OF POINTS
//MAKE A RESERVE SLIDE OUT WITH THE NON-POINT PLOTS



