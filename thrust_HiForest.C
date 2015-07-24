// Jennifer Coulter
// July 22th 2015
// Rutgers University, jennifer.coulter@cern.ch
//
// Test macro for plotting thrust, an event shape variable.
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

//plane class
class Plane{
public:
  TVector3 v1, v2;
  Plane(TVector3);
  
  //returns a projection onto the 2D plane 
  TVector3 Projection(TVector3 jaxis){
    //Find the projection of a jet onto this subspace
    Double_t scalar1 = jaxis.Dot(v1)/(v1.Dot(v1));
    Double_t scalar2 = (jaxis.Dot(v2)/v2.Dot(v2)); 
    TVector3 u1 = v1;   u1 = scalar1*u1;    
    TVector3 u2 = v2;   u2 = scalar2*u2;
    TVector3 proj = u1.operator+=(u2);
    return proj;
  }//end of projection
};
//plane class constructor
Plane::Plane(TVector3 nT){
  
  //Use TVector3 to find an orthogonal vector and a second vector orthogonal to the first and nT
  v1 = nT.Orthogonal();  v2 = nT.Cross(v1);

  //Normalize
  Double_t mag1 = v1.Mag();       Double_t mag2 = v2.Mag();
  v1(0) = v1(0)/mag1;    v1(1) = v1(1)/mag1;    v1(2) = v1(2)/mag1;
  v2(0) = v2(0)/mag2;    v2(1) = v2(1)/mag2;    v2(2) = v2(2)/mag2;	    
}//end plane constructor

//creates histograms in terms of Thrust vs. dN/dT
TH1F* DivideByBinWidth(TH1F * hist, const char * name){

  TH1F* h_return = new TH1F(name, "", hist->GetNbinsX(), 0,1);
  hist->Sumw2(); 
  //loops through all the bins
  for (int i=1;i<=hist->GetNbinsX();++i){
    Float_t bin = hist->GetBinWidth(i);
    Float_t val = hist->GetBinContent(i);
    Float_t valErr = h_return->GetBinError(i);
    val = val/bin;
    valErr= valErr/bin;
    h_return->SetBinError(i,valErr);
    h_return->SetBinContent(i, val); 
  }//end bin loop
  return h_return;
}//end rebin function

//Function to normalize a vector
TVector3 Norm(TVector3 v){
  if ( v(0) == 0 && v(1) == 0 && v(2) == 0) return v; 
  Double_t mag = TMath::Sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2)); 
  v(0) = v(0)/mag;    v(1) = v(1)/mag;   v(2) = v(2)/mag;
  return v; 
}//end normalize

//plot thrust
//void thrust_HiForest(Int_t startfile, Int_t endfile, Int_t jobNumber){
void thrust_HiForest(Int_t startfile = 10, Int_t endfile = 11, Int_t jobNumber = 1){

  TH1::SetDefaultSumw2();

  TStopwatch timer;
  timer.Start();
  
  bool debug = true;
  Float_t pT_cut = 30;
  Int_t radius = 3;
  jobNumber = 0; 
  
  //define trees and file
  //********* FIX FILE OPENING HERE
  TFile * file; 
  TFile * weight_file;
  TTree * t;
  TTree * hiEvt;
  TTree * hlt;
  TTree * skim;
  TTree * weight;
  TTree * thrust_tree;

  /*
  //set up the tree to be filled
  thrust_tree->Branch("pthatweight",&pthatweight,"pthatweight/D");
  thrust_tree->Branch("hiBin",&hiBin,"hiBin/I");
  thrust_tree->Branch("evt",&evnt,"evt/I");
  thrust_tree->Branch("lumi",&lumi,"lumi/I");
  thrust_tree->Branch("vz",&vz,"vz/F");
  */
  
  TH1F * h_thrust = new TH1F("thrust_unscaled", "", 50,0,1);
  TH1F * h_min = new TH1F("thrust_min", "", 50,0,1);
  TH1F * h_maj = new TH1F("thrust_maj", "", 50,0,1);
  TH1F * h_pT = new TH1F("pT", "", 100, 0, 120);
  TH1F * h_pTcut = new TH1F("pTcut", "", 100, 0, 120);
  TH1F * h_40 = new TH1F("thrust_40", "", 50,0,1);
  TH1F * h_60 = new TH1F("thrust_60", "", 50,0,1);
  TH1F * h_80 = new TH1F("thrust_80", "", 50,0,1);
  TH1F * h_nref = new TH1F("nref", "", 12, 0, 12);
  TH1F * h_jetCount = new TH1F("jetCount", "", 12, 0, 12);
  TH1F * h_eta = new TH1F("eta", "", 60, -2, 2);
  TH1F * h_phi = new TH1F("phi", "", 60, -3.15, 3.15);
  TH1F * h_weight = new TH1F("weighting", "", 500, 0, 500);

  //Tree variables
  Double_t px[1000];     Float_t pt[1000];    Int_t jt80;   Int_t jt80_pre;
  Double_t py[1000];     Float_t eta[1000];   Int_t jt40;   Int_t jt60_pre;
  Double_t pz[1000];     Float_t phi[1000];   Int_t jt60;   Int_t jt40_pre;
  
  Int_t nref;
  Float_t vz;
  Int_t halo;
  Int_t noise;
  Double_t pThat_weight;
  Float_t pThat; 

  Float_t dot = 0;
  Double_t mag = 0;
  Double_t thrust_temp = 0;
  Double_t thrust_max = 0;
  Double_t dot_maj = 0;
  Double_t dot_min = 0;
  Double_t min_temp = 0;
  Double_t maj_temp = 0;
  Double_t thrust_maj_max =0;
  Double_t thrust_min_max = 0;
  TVector3 max_thrust_axis;
  TVector3 p3Norm;
  Float_t max_eta = 0;   Float_t temp_eta = 0;  
  Float_t max_phi = 0;   Float_t temp_phi = 0;  
  Int_t jetCount = 0; //in order to sum up the number of jets per event
  Int_t eventCount = 0;//check to see how many of the events in each file are actually being used

  //define instream
  string input_file = "pp_MC_HiForest.txt"; 
  ifstream count(input_file.c_str(), ifstream::in);
  Int_t fileCount = 0;
  string * filename = new string[135];

  //count up the total number of files and save their names to an array of filenames for initialization purposes
  
  string line;
  while(getline(count, line)){
    filename[fileCount] = line;
    if (debug) cout << filename[fileCount] << endl; 
    fileCount++;
  }
  count.close();

  // For every file in file list, process trees
  for(int ifile = startfile; ifile < endfile; ifile++){
    
    string s = filename[ifile];
    string w = Form("weights/weights_pp_%d.root", ifile+1); 
    file = TFile::Open(s.c_str());
    weight_file = TFile::Open(w.c_str());     

    if (debug) cout << "\n **** =========================== New File ================================= **** \n ";
    cout << "File Name: " << filename[ifile] << endl; 
    cout << "File Number: " << endfile-ifile << "/" << endfile-startfile << endl;
    cout << "Weight File: " << w << endl;

    //define trees and file
    t = (TTree*)file->Get(Form("ak%dPFJetAnalyzer/t", radius));
    hiEvt = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
    hlt = (TTree*)file->Get("hltanalysis/HltTree");
    skim = (TTree*)file->Get("skimanalysis/HlTree");
    weight = (TTree*)weight_file->Get("weights");
    
    //Set branches of the tree 
    t->SetBranchAddress("jtpt", &pt);
    t->SetBranchAddress("jteta", &eta);
    t->SetBranchAddress("jtphi", &phi);
    t->SetBranchAddress("nref", &nref);
    t->SetBranchAddress("pthat", &pThat);

    hiEvt->SetBranchAddress("vz", &vz);

    //skim->SetBranchAddress("pHBHENoiseFilter", &noise);
    //skim->SetBranchAddress("pPAcollisionEventSelectionPA",&halo);
  
    hlt->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jt80);
    hlt->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jt60);
    hlt->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jt40);
    hlt->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jt80_pre);
    hlt->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jt60_pre);
    hlt->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jt40_pre);

    weight->SetBranchAddress("pthatweight", &pThat_weight);

    t->AddFriend(hiEvt);
    t->AddFriend(skim);
    t->AddFriend(hlt);
    t->AddFriend(weight);
    
    Long64_t nentries = t->GetEntries();
    nentries = 1;

    cout << "Events in File: " << nentries << endl;
    eventCount = 0;
    
    //event loop
    for(Long64_t nentry = 0; nentry<nentries; ++nentry){

      if(nentry%10000 == 0) cout << nentry << "%" << endl;
      
      t->GetEntry(nentry);
      jetCount = 0;
      bool select = false;

      //fill pThat spectra plot
      h_weight->Fill(pThat, pThat_weight); 

      //make selection cuts
      //if((TMath::Abs(vz) > 15)||(noise==0)||(halo==0)) {continue;}
      if((TMath::Abs(vz) > 15)) {continue;}
    
      //intial pT cut
      for(int k = 0; k < nref; k++){
	if(pt[k] >  30){
	  if((TMath::Abs(eta[k])<2)) select = true; 
	}
      }
      if(!select) continue;
      if(debug) cout<< " \n ******* New Event ******** " << endl;

      //reset maximum values
      eventCount++;
      thrust_max = 0; 

      //Part 1: Runs through all the jets in an event, checking them to see if they are the axis that maximizes thrust
      //max axis finding loop
      for(Long64_t naxis = 0; naxis < nref; ++naxis){

	//Cut checks
	if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[naxis]<<"\n \t eta = "<<eta[naxis]<<"\n \t phi = "<<phi[naxis]<<endl;
	h_pT->Fill(pt[naxis]);
	if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)) {
	  continue;
	}
	h_pTcut->Fill(pt[naxis]);
      
	jetCount++; 
      
	if(debug) cout<< " \n --------- New Test Axis (Thrust)--------- " << endl; 

	//reset values for this particular event
	thrust_temp = 0;   maj_temp = 0;   min_temp = 0;
      
	px[naxis] = pt[naxis]*TMath::Cos(phi[naxis]);
	py[naxis] = pt[naxis]*TMath::Sin(phi[naxis]);
	pz[naxis] = pt[naxis]*TMath::SinH(eta[naxis]);
      
	//calculates axis for this particular jet
	TVector3 nT (px[naxis], py[naxis], pz[naxis]);
	nT = Norm(nT);

	if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[naxis]<<"\n \t eta = "<<eta[naxis]<<"\n \t phi = "<<phi[naxis]<<endl;
	temp_phi = phi[naxis];   temp_eta = eta[naxis]; 
      
	//resets for next jet loop
	dot = 0;   mag = 0;

	//Part 2: Loops through all the jets to find the thrust value for the chosen axis 
	//jet loop
	for(Long64_t njet = 0; njet < nref; ++njet){

	  if((pt[njet] < pT_cut)||(TMath::Abs(eta[njet]) > 2)){ continue;}
	
	  if(debug) cout<< " \n --------- New Jet (Thrust)--------- " << endl; 
	  if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[njet]<<"\n \t eta = "<<eta[njet]<<"\n \t phi = "<<phi[njet]<<endl;
	
	  //calculate px, py, pz
	  px[njet] = pt[njet]*TMath::Cos(phi[njet]);
	  py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
	  pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	
	  //define momentum three vector
	  TVector3 p3 (px[njet], py[njet], pz[njet]);
	  //TVector3 p3Norm = Norm(p3);
	
	  if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	
	  //dots the two vectors for Thrust, Tmin and Tmaj
	  dot += TMath::Abs(p3.Dot(nT)); 
	  if(debug) cout<<"dot sum = " << dot << endl;
	
	  //sum the total p from the individual p magnitudes
	  mag += TMath::Abs(p3.Mag());
	  if(debug) cout<<"mag sum = " << mag << endl;
	
	}//end jet loop
      
	//calculate the thrust
	thrust_temp = ((dot)/mag);

	//Compare to see if this axis is a new maximum 
	if(debug) cout<< "\ntemp thrust = " << thrust_temp << endl; 
      
	if(thrust_temp>thrust_max){
	  thrust_max = thrust_temp;
	  max_thrust_axis = nT;
	  max_eta = temp_eta;
	  max_phi = temp_phi; 
	}
      
	if(debug) cout<< "max thrust = " << thrust_max << endl;
      
      }//end axis loop

      //Part 3: Begin code to select the Thrust Major and Minor axes
 
      //define the plane perpendicular to this axis in order to calculate Tmaj and Tmin
      Plane* perp = new Plane(max_thrust_axis);

      //reset maximum values for new axis test
      thrust_maj_max = 0;   thrust_min_max = 0; 
    
      //Thrust maj axis loop
      for(Long64_t naxis = 0; naxis < nref; ++naxis){

	if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)){ continue;}
	if(debug) cout<< " \n --------- New Test Axis (Min/Maj)--------- " << endl; 
      
	//define the jet axis for this iteration
	//calculate px, py, pz
	px[naxis] = pt[naxis]*TMath::Cos(phi[naxis]);
	py[naxis] = pt[naxis]*TMath::Sin(phi[naxis]); 
	pz[naxis] = pt[naxis]*TMath::SinH(eta[naxis]);
      
	//define momentum three vector
	TVector3 p3 (px[naxis], py[naxis], pz[naxis]);
	TVector3 p3Norm = Norm(p3);
            
	//define maj_axis and min_axis 
	TVector3 maj_axis = perp->Projection(p3Norm);
	maj_axis = Norm(maj_axis);
	TVector3 min_axis = max_thrust_axis.Cross(maj_axis);
	min_axis = Norm(min_axis); 
      
	if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	if(debug) cout<<"Maj Axis = {" << maj_axis(0) << ", " << maj_axis(1) << ", " << maj_axis(2)<< "}" << endl;
	if(debug) cout<<"Min Axis = {" << min_axis(0) << ", " << min_axis(1) << ", " << min_axis(2)<< "}\n" << endl;

	//reset for new axis test
	dot_maj = 0;   dot_min = 0;   mag = 0;
      
	//Part 4: Test the axis defined by the above loop to determine if this axis is the maximum
	//jet loop
	for(Long64_t njet = 0; njet < nref; ++njet){

	  //make a ptcut
	  if((pt[njet] < pT_cut)||(TMath::Abs(eta[njet]) > 2)){ continue;}
	
	  if(debug) cout<< " \n --------- New Jet (Maj/Min)--------- " << endl; 
	  if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[njet]<<"\n \t eta = "<<eta[njet]<<"\n \t phi = "<<phi[njet]<<endl;
	
	  //calculate px, py, pz
	  px[njet] = pt[njet]*TMath::Cos(phi[njet]);
	  py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
	  pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	
	  //define momentum three vector
	  TVector3 p3 (px[njet], py[njet], pz[njet]);
	  TVector3 p3Norm = Norm(p3);
	
	  if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	
	  //dots the two vectors for Tmin and Tmaj
	  dot_maj += TMath::Abs(p3.Dot(maj_axis)); 
	  dot_min += TMath::Abs(p3.Dot(min_axis));
	  if(debug) cout<<"dot maj sum = " << dot_maj << endl;
	  if(debug) cout<<"dot min sum = " << dot_min << endl;
	
	  //sum the total p from the individual p magnitudes
	  mag += TMath::Abs(p3.Mag());
	  if(debug) cout<<"mag sum = " << mag << endl;
	
	}//end jet loop

	//calculate the thrust major and minor for this axis
	maj_temp = dot_maj/mag;
	min_temp = dot_min/mag;
      
	//test to to see if this particular Tmaj and Tmin are the new maxima
	if(maj_temp>thrust_maj_max){
	  thrust_maj_max = maj_temp;  
	  thrust_min_max = min_temp; 
	  if(debug) cout << "thrust major max = "<< thrust_maj_max<< endl; 
	}   
      }//end of major/minor axis loop

      //fill all the maximum values before finishing
      if(jetCount > 1){
	
	h_thrust->Fill(thrust_max, pThat_weight);
	h_eta->Fill(max_eta, pThat_weight);
	h_phi->Fill(max_phi, pThat_weight);
	h_min->Fill(thrust_min_max, pThat_weight);
	h_maj->Fill(thrust_maj_max, pThat_weight);
	h_nref->Fill(nref);
	h_jetCount->Fill(jetCount);

	if(debug) {
	  if (thrust_max < 0.5)       cout << "FLAG_thrust1: " << thrust_max <<  " , " << jetCount << endl; 
	  if (thrust_maj_max > 0.5)   cout << "FLAG_maj: " << thrust_maj_max <<  " , " << jetCount << endl; 
	  if (thrust_min_max > 0.5)   cout << "FLAG_min: " << thrust_min_max <<  " , " << jetCount << endl;
	}
	
	if(jt80)    h_80->Fill(thrust_max,jt80_pre * pThat_weight);
	if(jt60)    h_60->Fill(thrust_max,jt60_pre * pThat_weight);
	if(jt40)    h_40->Fill(thrust_max,jt40_pre * pThat_weight);
      }
    }//end of event loop

    gROOT->GetListOfFiles()->Remove(file);
    gROOT->GetListOfFiles()->Remove(weight);
    
    cout << "Events Selected: " << eventCount << endl;
    cout << "File Finished" << endl; 
    
  }//end file loop

  //scale the histograms
  Double_t integral;
  
  integral = h_thrust->Integral();
  h_thrust->Scale(integral);

  integral = h_maj->Integral(); 
  h_maj->Scale(integral);

  integral = h_min->Integral();
  h_min->Scale(integral);

  integral = h_80->Integral();
  h_80->Scale(integral);

  integral = h_60->Integral();
  h_60->Scale(integral);

  integral = h_40->Integral();
  h_40->Scale(integral); 
 
  //Create the plot for Thrust vs. dN/dT
  //define histograms
  TH1F * h_T = DivideByBinWidth(h_thrust, "thrust_scaled");
  TH1F * h_Tmaj = DivideByBinWidth(h_maj, "thrust_maj_scaled");
  TH1F * h_Tmin = DivideByBinWidth(h_min, "thrust_min_scaled");
  
  Float_t divEntries = 1./(h_thrust->GetBinWidth(1));
  
  h_T->Scale(divEntries);
  h_Tmaj->Scale(divEntries);
  h_Tmin->Scale(divEntries);

  h_40 = DivideByBinWidth(h_40, "thrust_40_new");   h_40->Scale(divEntries);
  h_60 = DivideByBinWidth(h_60, "thrust_60_new");   h_60->Scale(divEntries);
  h_80 = DivideByBinWidth(h_80, "thrust_80_new");   h_80->Scale(divEntries);

  TFile * save_File = new TFile(Form("test_pp_thrust_%d.root", jobNumber),"RECREATE");
  
  h_T->SetDirectory(save_File);
  h_Tmaj->SetDirectory(save_File); 
  h_Tmin->SetDirectory(save_File); 
  h_pT->SetDirectory(save_File); 
  h_pTcut->SetDirectory(save_File); 
  h_40->SetDirectory(save_File); 
  h_60->SetDirectory(save_File); 
  h_80->SetDirectory(save_File);
  h_nref->SetDirectory(save_File); 
  h_jetCount->SetDirectory(save_File); 
  h_eta->SetDirectory(save_File); 
  h_phi->SetDirectory(save_File);
  h_weight->SetDirectory(save_File); 
  
  h_T->Write();
  h_Tmaj->Write();
  h_Tmin->Write();
  h_pT->Write();
  h_pTcut->Write();
  h_40->Write();
  h_60->Write();
  h_80->Write();
  h_nref->Write();
  h_jetCount->Write();
  h_eta->Write();
  h_phi->Write();
  h_weight->Write(); 
  
  save_File->Write();
  save_File->Close();

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
    
  }//end of plot thrust
