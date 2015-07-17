// Jennifer Coulter
// June 24th 2015
// Rutgers University, jennifer.coulter@cern.ch
//
//Macro to compute thrust for a large data set.
//Modified version of test_thrust.C

//July 16th -> added fix so that jettree and weight and event trees are on the same event number

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
    //cout<<"v1 = {" << v1(0) << ", " << v1(1) << ", " << v1(2)<< "}" << endl;
    //cout<<"v2 = {" << v2(0) << ", " << v2(1) << ", " << v2(2)<< "}" << endl;
    Double_t scalar1 = jaxis.Dot(v1)/(v1.Dot(v1));
    Double_t scalar2 = (jaxis.Dot(v2)/v2.Dot(v2)); 
    TVector3 u1 = v1;   u1 = scalar1*u1;    
    TVector3 u2 = v2;   u2 = scalar2*u2;
    TVector3 proj = u1.operator+=(u2);
    return proj;
  }//end of projection
  /*
  //returns a projection onto the 2D plane 
  TVector3 Projection(TVector3 jaxis){
    //Find the projection of a jet onto this subspace
    scalar1 = jaxis.Dot(v1)/(v1.Dot(v1));
    scalar2 = jaxis.Dot(v2)/(v2.Dot(v2)); 

    v1 = scalar1*v1;
    v2 = scalar2*v2;
    
    proj = v1.operator+=(v2);
    return proj;
  }//end of projection
  */
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
  //check for a zero vector
  if ( v(0) == 0 && v(1) == 0 && v(2) == 0) return v; 
  Double_t mag = TMath::Sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2)); 
  v(0) = v(0)/mag;    v(1) = v(1)/mag;   v(2) = v(2)/mag;
  return v; 
}//end normalize

//plot thrust
void thrust_data(Float_t pT_cut = 30, Int_t radius = 2, int startfile = 21, int endfile = 22){

  TH1::SetDefaultSumw2();
  
  TStopwatch timer;
  timer.Start();
  
  bool debug = false;

  Double_t px[1000];     Float_t  pt[1000];    Int_t jt80;   Int_t jt80_pre;
  Double_t py[1000];     Float_t eta[1000];   Int_t jt40;   Int_t jt60_pre;
  Double_t pz[1000];     Float_t phi[1000];   Int_t jt60;   Int_t jt40_pre;
  
  Int_t nref;
  Int_t isMultiMatch;
  Int_t isGoodEvt;
  Double_t pThat_weight;
  Float_t pThat; 
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
  TVector3 p3Norm;
  Float_t max_eta = 0;   Float_t temp_eta = 0;  
  Float_t max_phi = 0;   Float_t temp_phi = 0;
  Int_t jetTreeCount = 0; //a counter to make sure that the weight tree and jet tree are on the same 
  Int_t jetCount = 0; //in order to sum up the number of jets per event
  Int_t eventCount = 0;//check to see how many of the events in each file are actually being used

  //define instream
  string input_file = "pp_MC_ak2_ordered.txt"; 
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

  TTree * jet;
  TTree * event;
  TTree * weight;
  TTree * thrust_tree = new TTree("t_thrust", "Tree of thrust events for further analysis"); 
  TFile * file;
  TFile * weight_file;

  /*
  //set up the tree to be filled
  thrust_tree->Branch("pthatweight",&pthatweight,"pthatweight/D");
  thrust_tree->Branch("hiBin",&hiBin,"hiBin/I");
  thrust_tree->Branch("evt",&evnt,"evt/I");
  thrust_tree->Branch("lumi",&lumi,"lumi/I");
  thrust_tree->Branch("vz",&vz,"vz/F");
  */

  // For every file in file list, process trees
  for(int ifile = 0; ifile < 45; ifile++){
    
    string s = "root://eosuser.cern.ch://eos/user/j/jcoulter/Data/pp_MC/";
    string str = s.append(filename[ifile]);
    string w = Form("weights/weights_pp_%d.root", ifile+1); 
    file = TFile::Open(str.c_str());
    weight_file = TFile::Open(w.c_str());     

    if (debug) cout << "\n **** =========================== New File ================================= **** \n ";
    cout << "File Name: " << filename[ifile] << endl; 
    cout << "File Number: " << ifile << "/" << fileCount << endl;
    cout << "Weight File: " << w << endl;
    
    jet = (TTree*)file->Get(Form("ak%dJetAnalyzer/jetTree", radius));
    event = (TTree*)file->Get(Form("ak%dJetAnalyzer/evtTree", radius));
    weight = (TTree*)weight_file->Get("weights");
   
    //Set branches of the tree
    jet->SetBranchAddress("pfpt", &pt);
    jet->SetBranchAddress("pfeta", &eta);
    jet->SetBranchAddress("pfphi", &phi);
    jet->SetBranchAddress("npf", &nref);
    jet->SetBranchAddress("pthat", &pThat);

    Int_t run_J, lumi_J, evt_J;
    Float_t vz_J;
    jet->SetBranchAddress("run_value", &run_J);
    jet->SetBranchAddress("lumi_value", &lumi_J);
    jet->SetBranchAddress("evt_value", &evt_J);
    jet->SetBranchAddress("vz", &vz_J);

    //jet->SetBranchAddress("pPAcollisionEventSelectionPA",&halo);

    jet->SetBranchAddress("jet80",&jt80);   jet->SetBranchAddress("jet80_prescl",&jt80_pre);
    jet->SetBranchAddress("jet60",&jt60);   jet->SetBranchAddress("jet60_prescl",&jt60_pre);
    jet->SetBranchAddress("jet40",&jt40);   jet->SetBranchAddress("jet40_prescl",&jt40_pre);

    //branches used to check if the events have the same value
    Int_t run_E, lumi_E, evt_E;
    Float_t vz_E; 
    event->SetBranchAddress("run_value", &run_E);
    event->SetBranchAddress("lumi_value", &lumi_E);
    event->SetBranchAddress("evt_value", &evt_E);
    event->SetBranchAddress("vz", &vz_E);
    event->SetBranchAddress("isGoodEvt", &isGoodEvt);
   
    weight->SetBranchAddress("pthatweight", &pThat_weight);

    jet->AddFriend(event);
    jet->AddFriend(weight);
    
    Long64_t nentries = event->GetEntries();
    //if(debug) nentries = 1000;
    nentries = 100000; 
    
    cout << "Events in File: " << nentries << endl;
    eventCount = 0;
    jetTreeCount = 0;
    
    //if(event->GetEntries() != weight->GetEntries()) cout<<"Event tree and weight tree have different entries"<<endl;
    //if(jet->GetEntries() != weight->GetEntries()) cout<<"Jet tree and weight tree have different entries"<<endl;
    // if(jet->GetEntries() != event->GetEntries()) cout<<"Jet tree and event tree have different entries"<<endl;

    //cout<<"jet Tree events = "<<jet->GetEntries()<<endl;
    //cout<<"event Tree events = "<<event->GetEntries()<<endl;
    //cout<<"weight Tree events = "<<weight->GetEntries()<<endl;
    
    //event loop
    for(Long64_t nentry = 0; nentry<nentries; ++nentry){

      jet->GetEvent(jetTreeCount);
      event->GetEvent(nentry);
      jetCount = 0; 
      bool select = false;

      if(nentry%10000 == 0) cout << nentry/nentries << "%" << endl;
      
      //make selection cuts -> if event tree is not a good event, trash it
      //if((TMath::Abs(vz_E) > 15)||(!isGoodEvt)) {continue;}
      if(!isGoodEvt) {continue;}

      if (debug) cout << "event: " << "vz: " << vz_E << " lumi: " << lumi_E << " event: " << evt_E << endl;
      if (debug) cout << "jet: " << "vz: " << vz_J << " lumi: " << lumi_J << " event: " << evt_J << endl;

      //check if jet and weight are on the same event, if they are the same, increment the jetTreeCount variable to keep weight and jet trees in sync
      if( (vz_J != vz_E) && (lumi_J != lumi_E) && (evt_J != evt_E)) {  cout << "RED FLAG!!!! event not same, JetTreeCount: " << jetTreeCount << " nentry: " << nentry << endl;  continue;}
      jetTreeCount++; 
      
      //fill pThat spectra plot
      h_weight->Fill(pThat, pThat_weight); 
      
      //intial pT cut
      for(int k = 0; k < nref; k++){
	if(pt[k] >  pT_cut){
	  if((TMath::Abs(eta[k])<2)) select = true; 
	}
      }
      if(!select) continue;
      if(debug) cout << " \n ******* New Event ******** " << endl;
      if(debug) cout << "nref: " << nref << endl;
   
      //reset maximum values
      eventCount++;
      thrust_max = 0; 
      //Part 1: Runs through all the jets in an event, checking them to see if they are the axis that maximizes thrust
      //max axis finding loop
      for(Long64_t naxis = 0; naxis < nref; ++naxis){

	//Cut checks
	h_pT->Fill(pt[naxis], pThat_weight);

	if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)) {continue;}	
	h_pTcut->Fill(pt[naxis], pThat_weight);
			
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
	  
	  if((pt[njet] < pT_cut)||(TMath::Abs(eta[njet]) > 2)){continue;}
	  
	  if(debug) cout<< " \n --------- New Jet (Thrust)--------- " << endl; 
	  if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[njet]<<"\n \t eta = "<<eta[njet]<<"\n \t phi = "<<phi[njet]<<endl;
	  
	  //calculate px, py, pz
	  px[njet] = pt[njet]*TMath::Cos(phi[njet]);
	  py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
	  pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	  
	  //define momentum three vector
	  TVector3 p3jet (px[njet], py[njet], pz[njet]);
	  //p3Norm = Norm(p3);
	  
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
	if(debug) cout << "max thrust = " << thrust_max << endl;

      }//end axis loop
 
      //Part 3: Begin code to select the Thrust Major and Minor axes

      //define the plane perpendicular to this axis in order to calculate Tmaj and Tmin
      Plane * perp = new Plane(max_thrust_axis);
      
      //reset maximum values for new axis test
      thrust_maj_max = 0;   thrust_min_max = 0; 
      
      //Thrust maj axis loop
      for(Long64_t naxis = 0; naxis < nref; ++naxis){
	
	if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)){continue;}
	if(debug) cout<< " \n --------- New Test Axis (Min/Maj)--------- " << endl; 
	
	//define the jet axis for this iteration
	//calculate px, py, pz
	px[naxis] = pt[naxis]*TMath::Cos(phi[naxis]);
	py[naxis] = pt[naxis]*TMath::Sin(phi[naxis]); 
	pz[naxis] = pt[naxis]*TMath::SinH(eta[naxis]);
	
	//define momentum three vector
	TVector3 p3 (px[naxis], py[naxis], pz[naxis]);
	
	//define maj_axis and min_axis 
	TVector3 maj_axis = perp->Projection(Norm(p3));
	if(debug) cout<<"Maj Axis = {" << maj_axis(0) << ", " << maj_axis(1) << ", " << maj_axis(2)<< "}" << endl;
	maj_axis = Norm(maj_axis);
	TVector3 min_axis = max_thrust_axis.Cross(maj_axis);
	min_axis = Norm(min_axis); 
	
	if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	if(debug) cout<<"Thr Axis = {" << max_thrust_axis(0) << ", " << max_thrust_axis(1) << ", " << max_thrust_axis(2) <<  "}" << endl;
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
	  if(debug) cout<< "Jet Variables: " << "\n \t pT = " << pt[njet] << "\n \t eta = " << eta[njet] << "\n \t phi = " << phi[njet] << endl;
	  
	  //calculate px, py, pz
	  px[njet] = pt[njet]*TMath::Cos(phi[njet]);
	  py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
	  pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	  
	  //define momentum three vector
	  TVector3 p3jet (px[njet], py[njet], pz[njet]);
	  //p3Norm = Norm(p3);
	  
	  if(debug) cout << "Jet Axis = {" << p3jet(0) << ", " << p3jet(1) << ", " << p3jet(2)<< "}" << endl;
	  
	  //dots the two vectors for Tmin and Tmaj
	  if(debug) cout<<"Maj Axis = {" << maj_axis(0) << ", " << maj_axis(1) << ", " << maj_axis(2)<< "}" << endl;
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
	  if(debug) cout << "thrust minor max = "<< thrust_min_max<< endl; 
	}   

      }//end of major/minor axis loop

      //fill all the maximum values before finishing, unless there is only one selected jet
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
	
	if(jt80)   h_80->Fill(thrust_max,jt80_pre*pThat_weight);
	if(jt60)   h_60->Fill(thrust_max,jt60_pre*pThat_weight);
	if(jt40)   h_40->Fill(thrust_max,jt40_pre*pThat_weight);
      } 
      
    }//end of event loop

    file->Close();
    //TROOT::gROOT->GetListOfFiles()->Remove(file);
    //ROOT.gROOT.GetListOfFiles().Remove(weight);
    weight_file->Close();

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
  h_weight->SetDirectory(save_File); 

  h_T->Write();
  h_Tmaj->Write();
  h_Tmin->Write();
  h_thrust->Write();
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


//for each file create an event level output tree with seven branches
//vs run lumi event thrust thrust min thrust maj
//add tree to an output file
//can look at thrust from .8 to 1 for a certain type of event like dijet
