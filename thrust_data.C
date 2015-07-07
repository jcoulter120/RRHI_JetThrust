// Jennifer Coulter
// June 24th 2015
// Rutgers University, jennifer.coulter@cern.ch
//
//Macro to compute thrust for a large data set.
//Modified version of test_thrust.C

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

//plane class
class Plane{
public:
  TVector3 v1, v2, u1, u2, proj;
  Double_t scalar1, scalar2, mag1, mag2; 
  Plane(TVector3);
  
  //returns a projection onto the 2D plane 
  TVector3 Projection(TVector3 jaxis){
    //Find the projection of a jet onto this subspace
    scalar1 = jaxis.Dot(v1)/(v1.Dot(v1));
    scalar2 = jaxis.Dot(v2)/(v2.Dot(v2)); 
    u1 = v1;   u1 = scalar1*u1;    
    u2 = v2;   u2 = scalar2*u2;
    proj = u1.operator+=(u2);
    return proj;
  }//end of projection
};
//plane class constructor
Plane::Plane(TVector3 nT){
  
  //Use TVector3 to find an orthogonal vector and a second vector orthogonal to the first and nT
  v1 = nT.Orthogonal();  v2 = nT.Cross(v1);

  //Normalize
  mag1 = v1.Mag();       mag2 = v2.Mag();
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
    Float_t valErr =hist->GetBinError(i);
    val = val/bin;
    valErr= valErr/bin;
    h_return->SetBinError(i,valErr);
    h_return->SetBinContent(i, val); 
  }//end bin loop
  return h_return;
}//end rebin function

//Function to normalize a vector
TVector3 Norm(TVector3 v){
  Double_t mag = TMath::Sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2)); 
  v(0) = v(0)/mag;    v(1) = v(1)/mag;   v(2) = v(2)/mag;
  return v; 
}//end normalize

//plot thrust
void thrust_data(Float_t pT_cut = 15, Int_t radius = 2, int startfile = 21, int endfile = 22){

  TH1::SetDefaultSumw2();
  bool debug = false;

  //TFile * save_File = new TFile("test_pp_thrust.root","RECREATE");

  Double_t px[10000];     Float_t  pt[10000];    Int_t jt80;   Int_t jt80_pre;
  Double_t py[10000];     Float_t eta[10000];   Int_t jt40;   Int_t jt60_pre;
  Double_t pz[10000];     Float_t phi[10000];   Int_t jt60;   Int_t jt40_pre;
  
  Int_t nref;
  Float_t vz;
  Int_t isMultiMatch;
  Int_t isGoodEvt;
  //Int_t halo;
  //Int_t noise;

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
  Float_t max_eta = 0;   Float_t temp_eta = 0;  
  Float_t max_phi = 0;   Float_t temp_phi = 0;  
  Int_t jetCount = 0; //in order to sum up the number of jets per event

  //define instream
  string input_file = "pp_MC_ak2.txt"; 
  ifstream count(input_file.c_str(), ifstream::in);
  Int_t total_fileCount = 0;
  Int_t filesUsed = 0; //this is the number of files used in the current set of histograms
  string * filename = new string[135];

  //count up the total number of files and save their names to an array of filenames for initialization purposes

  string line;
  while(getline(count, line)){
    filename[total_fileCount] = line;
    if (debug) cout << filename[total_fileCount] << endl; 
    total_fileCount++;
  }
  count.close();
  
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

  TTree * jet;
  TTree * event;
  TFile * file;

  //loop to initialize all the files
  //ifstream in(input_file.c_str(), ifstream::in);

  // For every file in file list, process trees
  //for(int ifile = 0; ifile < 45; ifile++){
  for(int ifile = 0; ifile < 10; ifile++){

    cout << "made it here" << endl; 
    
    string s = "root://eosuser.cern.ch://eos/user/j/jcoulter/Data/pp_MC/";
    string str = s.append(filename[ifile]);
    file = TFile::Open(str.c_str()); 
  //   //cout << filename[ifile] << endl;
  //   //cout << str << endl; 
   
    if (debug) cout << "\n **** =========================== New File ================================= **** \n ";
    cout << "File Name: " << filename[ifile] << endl; 
    cout << "File Number: " << ifile << "/" << total_fileCount << endl;
    if (debug) cout << "actual file name = " << str.c_str() <<endl;
    
    jet = (TTree*)file->Get(Form("ak%dJetAnalyzer/jetTree", radius));
    event = (TTree*)file->Get(Form("ak%dJetAnalyzer/evtTree", radius));
    
    //Set branches of the tree
    jet->SetBranchAddress("pfpt", &pt);
    jet->SetBranchAddress("pfeta", &eta);
    jet->SetBranchAddress("pfphi", &phi);
    jet->SetBranchAddress("npf", &nref);

    jet->AddFriend(event);
    jet->SetBranchAddress("vz", &vz);
    jet->SetBranchAddress("isMultiMatch", &isMultiMatch);
    jet->SetBranchAddress("isGoodEvt", &isGoodEvt);
  
    //t->SetBranchAddress("pPAcollisionEventSelectionPA",&halo); 
    //t->SetBranchAddress("pHBHENoiseFilter", &noise);
    
    jet->SetBranchAddress("jet80",&jt80);   jet->SetBranchAddress("jet80_prescl",&jt80_pre);
    jet->SetBranchAddress("jet60",&jt60);   jet->SetBranchAddress("jet60_prescl",&jt60_pre);
    jet->SetBranchAddress("jet40",&jt40);   jet->SetBranchAddress("jet40_prescl",&jt40_pre);

    Long64_t nentries = jet->GetEntries();
    cout << "Events in File: " << nentries << endl; 
    if(debug) nentries = 100;
  
    //event loop
    for(Long64_t nentry = 0; nentry<10; ++nentry){
      
      jet->GetEntry(nentry);
      jetCount = 0;
      bool select = false;
      
      //make selection cuts
      //if((TMath::Abs(vz) > 15)||(noise==0)||(halo==0)||(!isGoodEvt) {continue;}
      if((TMath::Abs(vz) > 15)||(!isGoodEvt)) {continue;}
      
      //intial pT cut
      for(int k = 0; k < nref; k++){
	if(pt[k] >  30){
	  if((TMath::Abs(eta[k])<2)) select = true; 
	}
      }
      if(!select) continue;
      if(debug) cout<< " \n ******* New Event ******** " << endl;
      
      //reset maximum values
      thrust_max = 0; 
      
      //Part 1: Runs through all the jets in an event, checking them to see if they are the axis that maximizes thrust
      //max axis finding loop
      for(Long64_t naxis = 0; naxis < nref; ++naxis){

	//Cut checks
	if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[naxis]<<"\n \t eta = "<<eta[naxis]<<"\n \t phi = "<<phi[naxis]<<endl;
	h_pT->Fill(pt[naxis]);
	if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)||(isMultiMatch)) {continue;}
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
	  
	  if((pt[njet] < pT_cut)||(TMath::Abs(eta[njet]) > 2)||(isMultiMatch)){continue;}
	  
	  if(debug) cout<< " \n --------- New Jet (Thrust)--------- " << endl; 
	  if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[njet]<<"\n \t eta = "<<eta[njet]<<"\n \t phi = "<<phi[njet]<<endl;
	  
	  //calculate px, py, pz
	  px[njet] = pt[njet]*TMath::Cos(phi[njet]);
	  py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
	  pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	  
	  //define momentum three vector
	  TVector3 p3jet (px[njet], py[njet], pz[njet]);
	  //TVector3 p3Norm = Norm(p3);
	  
	  if(debug) cout<<"Jet Axis = {" << p3jet(0) << ", " << p3jet(1) << ", " << p3jet(2)<< "}" << endl;
	
	  //dots the two vectors for Thrust, Tmin and Tmaj
	  dot += TMath::Abs(p3jet.Dot(nT)); 
	  if(debug) cout<<"dot sum = " << dot << endl;
	  
	  //sum the total p from the individual p magnitudes
	  mag += TMath::Abs(p3jet.Mag());
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
	
	if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)||(isMultiMatch)){continue;}
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
	TVector3 min_axis = max_thrust_axis.Cross(maj_axis);
	
	if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	if(debug) cout<<"Maj Axis = {" << maj_axis(0) << ", " << maj_axis(1) << ", " << maj_axis(2)<< "}" << endl;
	if(debug) cout<<"Min Axis = {" << min_axis(0) << ", " << min_axis(1) << ", " << min_axis(2)<< "}\n" << endl;
	
	//reset for new axis test
	dot_maj = 0;   dot_min = 0;   mag = 0;
	
	//Part 4: Test the axis defined by the above loop to determine if this axis is the maximum
	//jet loop
	for(Long64_t njet = 0; njet < nref; ++njet){
	  
	  //make a ptcut
	  if((pt[njet] < pT_cut)||(TMath::Abs(eta[njet]) > 2)||(isMultiMatch)){ continue;}
	  
	  if(debug) cout<< " \n --------- New Jet (Maj/Min)--------- " << endl; 
	  if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[njet]<<"\n \t eta = "<<eta[njet]<<"\n \t phi = "<<phi[njet]<<endl;
	  
	  //calculate px, py, pz
	  px[njet] = pt[njet]*TMath::Cos(phi[njet]);
	  py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
	  pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	  
	  //define momentum three vector
	  TVector3 p3jet (px[njet], py[njet], pz[njet]);
	  //TVector3 p3Norm = Norm(p3);
	  
	  if(debug) cout<<"Jet Axis = {" << p3jet(0) << ", " << p3jet(1) << ", " << p3jet(2)<< "}" << endl;
	  
	  //dots the two vectors for Tmin and Tmaj
	  dot_maj += TMath::Abs(p3jet.Dot(maj_axis)); 
	  dot_min += TMath::Abs(p3jet.Dot(min_axis));
	  if(debug) cout<<"dot maj sum = " << dot_maj << endl;
	  if(debug) cout<<"dot min sum = " << dot_min << endl;
	  
	  //sum the total p from the individual p magnitudes
	  mag += TMath::Abs(p3jet.Mag());
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
      
      //fill all the maximum values before finishing, unless there is only one selected jet
      if(jetCount != 1){
	
	h_thrust->Fill(thrust_max);
	h_eta->Fill(max_eta);
	h_phi->Fill(max_phi);
	h_min->Fill(thrust_min_max);
	h_maj->Fill(thrust_maj_max);
	h_nref->Fill(nref);
	h_jetCount->Fill(jetCount);

	if(jt80)   h_80->Fill(thrust_max,jt80_pre);
	if(jt60)   h_60->Fill(thrust_max,jt60_pre);
	if(jt40)   h_40->Fill(thrust_max,jt40_pre);
      }    
    }//end of event loop
    
    h_thrust->Scale(nentries);
    h_maj->Scale(nentries);
    h_min->Scale(nentries);
    h_80->Scale(nentries);
    h_60->Scale(nentries);
    h_40->Scale(nentries); 

    file->Close();
    filesUsed++;
    cout << "file finished" << endl; 
    
  }//end file loop

  //in.close(); 
  //Create the plot for Thrust vs. dN/dT
  //define histograms
  TH1F * h_T = DivideByBinWidth(h_thrust,"thrust_scaled");
  TH1F * h_Tmaj = DivideByBinWidth(h_maj, "thrust_maj_scaled");
  TH1F * h_Tmin = DivideByBinWidth(h_min, "thrust_min_scaled");
  
  Float_t divEntries = 1./(h_thrust->GetBinWidth(1));
  
  h_T->Scale(divEntries);
  h_Tmaj->Scale(divEntries);
  h_Tmin->Scale(divEntries);
  
  h_40 = DivideByBinWidth(h_40, "thrust_40_new");   h_40->Scale(divEntries);
  h_60 = DivideByBinWidth(h_60, "thrust_60_new");   h_60->Scale(divEntries);
  h_80 = DivideByBinWidth(h_80, "thrust_80_new");   h_80->Scale(divEntries);

  TFile * save_File = new TFile("test_pp_thrust.root","RECREATE");
  
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
  
  save_File->Close();
  
  cout<<"Macro End"<<endl;
}//end of plot thrust
