// Jennifer Coulter
// August 13th 2015
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
  TVector3 v1, v2, proj, u1, u2;
  Double_t scalar1, scalar2, mag1, mag2; 
  Plane(TVector3);

  //returns a projection onto the 2D plane 
  TVector3 Projection(TVector3 jaxis){
    //Find the projection of a jet onto this subspace
    if(v1.Mag() == 0) { scalar1 = 0; }   else { scalar1 = jaxis.Dot(v1)/(v1.Dot(v1)); } 
    if(v2.Mag() == 0) { scalar2 = 0; }   else { scalar2 = jaxis.Dot(v2)/(v2.Dot(v2)); } 
    v1 = scalar1*v1;
    v2 = scalar2*v2;
    proj(0) = v1(0) + v2(0);
    proj(1) = v1(1) + v2(1);
    proj(2) = v1(2) + v2(2); 
    
    return proj;
  }//end of projection
};
//plane class constructor
Plane::Plane(TVector3 nT){
  
  //Use TVector3 to find an orthogonal vector and a second vector orthogonal to the first and nT
  v1 = nT.Orthogonal();  v2 = nT.Cross(v1);

  //Normalize, checking for 0 length axes
  if ((v1(0) == 0) && (v1(1) == 0) && (v1(2) == 0)){  v1(0) = 0;    v1(1) = 0;    v1(2) = 0; }
  else { mag1 = v1.Mag();   v1(0) = v1(0)/mag1;    v1(1) = v1(1)/mag1;    v1(2) = v1(2)/mag1; } 
  if ((v2(0) == 0) && (v2(1) == 0) && (v2(2) == 0)){  v2(0) = 0;    v2(1) = 0;    v2(2) = 0; } 
  else { mag2 = v2.Mag();   v2(0) = v2(0)/mag2;    v2(1) = v2(1)/mag2;    v2(2) = v2(2)/mag2; }
}//end plane constructor

//creates histograms in terms of Thrust vs. dN/dT
TH1F* DivideByBinWidth(TH1F * hist, const char * name){

  TH1F* h_return = new TH1F(name, "", hist->GetNbinsX(), 0,1);
  hist->Sumw2(); 
  //loops through all the bins
  for (int i=1;i<=hist->GetNbinsX();++i){
    Float_t bin = hist->GetBinWidth(i);
    Float_t val = hist->GetBinContent(i);
    Float_t valErr = hist->GetBinError(i);
    val = val/bin;
    valErr= valErr/bin;
    h_return->SetBinError(i,valErr);
    h_return->SetBinContent(i, val); 
  }//end bin loop
  return h_return;
}//end rebin function

//Function to normalize a vector
TVector3 Norm(TVector3 v){
  if ( (v(0) == 0) && (v(1) == 0) && (v(2) == 0)) return v; 
  Double_t mag = TMath::Sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2)); 
  v(0) = v(0)/mag;    v(1) = v(1)/mag;   v(2) = v(2)/mag;
  return v; 
}//end normalize

//plot thrust
//void thrust_HiForest(Int_t startfile, Int_t endfile, Int_t jobNumber){

void thrust_Pythia(Int_t startfile = 0,
		   Int_t endfile = 12,
		   Int_t jobNumber = 60,
		   int radius = 3,
		   float ptCutLead = 60.0,
		   float ptCut = 30.0,
		   float etaCut = 2.0){

  TH1::SetDefaultSumw2();

  TStopwatch timer;
  bool debug = false;
  //Float_t pT_cut = 30;
  //Int_t radius = 3;
  
  //define trees and file
  TFile * file; 
  TFile * save_File = new TFile(Form("pythia_thrust_%d.root", jobNumber),"RECREATE");
  TFile * bad_File = new TFile(Form("pythia_bad_file_%d.root", jobNumber),"RECREATE");

  TTree * t;
  TTree * hi;
  TTree * nt;

  TH1F * h_thrust = new TH1F("thrust_unscaled", "", 50,0,1);
  TH1F * h_min = new TH1F("thrust_min", "", 50,0,1);
  TH1F * h_maj = new TH1F("thrust_maj", "", 50,0,1);
  TH1F * h_pT = new TH1F("pT", "", 100, 0, 500);
  TH1F * h_pTcut = new TH1F("pTcut", "", 100, 0, 500);
  //TH1F * h_40 = new TH1F("thrust_40", "", 50,0,1);
  //TH1F * h_60 = new TH1F("thrust_60", "", 50,0,1);
  //TH1F * h_80 = new TH1F("thrust_80", "", 50,0,1);
  TH1F * h_nref = new TH1F("nref", "", 12, 0, 12);
  TH1F * h_jetCount = new TH1F("jetCount", "", 12, 0, 12);
  TH1F * h_eta = new TH1F("eta", "", 60, -2, 2);
  TH1F * h_phi = new TH1F("phi", "", 60, -3.15, 3.15);
  TH1F * h_weight = new TH1F("weighting", "", 1200, 0, 1200);
  TH1F * h_unweight = new TH1F("unweighted", "", 1200, 0, 1200);
  
  TH1F * h_TBad = new TH1F("thrust_bad", "", 50,0,1);
  TH1F * h_TminBad = new TH1F("thrust_min_bad", "", 50,0,1);
  TH1F * h_TmajBad = new TH1F("thrust_maj_bad", "", 50,0,1);
  TH1F * h_TBadpT = new TH1F("TBadpT", "", 100, 0, 120);
  TH1F * h_TmajBadpT = new TH1F("TmajBadpT", "", 100, 0, 120);
  TH1F * h_TminBadpT = new TH1F("TminBadpT", "", 100, 0, 120);
  TH1F * h_etaBad = new TH1F("etaBad", "", 60, -2, 2);
  TH1F * h_phiBad = new TH1F("phiBad", "", 60, -3.15, 3.15);
  TH1F * h_nrefBad = new TH1F("nrefBad", "", 12, 0, 12);
  TH1F * h_jetCountBad = new TH1F("jetCountBad", "", 12, 0, 12);
  TH1F * h_weightBad = new TH1F("weightingBad", "", 700, 0, 700);
  TH1F * h_pthatBad = new TH1F("pthatBad", "", 500, 0, 500);
  TH1F * h_fileNum = new TH1F("fileNum", "", 45, 0, 45);

  //Tree variables
  Float_t pt[1000];    Int_t jt80;   Int_t jt80_pre;
  Float_t eta[1000];   Int_t jt40;   Int_t jt60_pre;
  Float_t phi[1000];   Int_t jt60;   Int_t jt40_pre;
  
  Int_t nref;
  Float_t vz;
  Int_t lumi;
  Int_t evt = 0; //this is used as the event number out of all the events handled
  Int_t isBadThrust = 0;
  Int_t isGoodEvt = 0; 
  Int_t halo;
  Int_t noise;
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
  Int_t max_nref;
  Int_t jetCount = 0; //in order to sum up the number of jets per event
  Int_t eventCount = 0;//check to see how many of the events in each file are actually being used
  Int_t badCount = 0; 
  //Double_t weightQ;
  Double_t weightSTAT;

  //PTHAT WEIGHTING SETUP======================================================
  //The sums of the cross sections for each bin, as obtained from Jewel output log files, are divided by the number of events in each bin to get the weight for each event
  double pthatBinning[] = {15,30,50,80,120,170,220,280,330,400,460,540};
  //                0           1           2      3       4           5            6        7        8        9          10
  //double xsSUM[] = {4878812.4,4218696.2,4659866.4,4049616.0,2934218.1,1695768.4,1302251.3,718081.6,578537.1,291444.26,168207.79};
  //for 120k
  double xsSTAT[] = {0.193,0,0.0009228,8.835e-05,9.752e-06,1.182e-06,2.314e-07,3.795e-08,1.22e-08,2.202e-09,6.641e-10,1.483e-10};
  //for 1.2M

  /*
  //set up the tree to be filled
  TTree * thrust_tree = new TTree("Thrust_Values","");
  thrust_tree->Branch("pthatweight",&weightQ,"pthatweight/D");
  //thrust_tree->Branch("hiBin",&hiBin,"hiBin/I");
  thrust_tree->Branch("evt",&evt,"evt/I");
  //thrust_tree->Branch("run",&run,"run/I");
  //thrust_tree->Branch("lumi",&lumi,"lumi/I");
  thrust_tree->Branch("vz",&vz,"vz/F");
  //thrust_tree->Branch("isGoodEvt",&isGoodEvt,"isGoodEvt/I");
  thrust_tree->Branch("isBadThrust",&isBadThrust,"isBadThrust/I"); // 1 if thrust is bad, otherwise 0
  thrust_tree->Branch("thrust",&thrust_max,"thrust/D");
  thrust_tree->Branch("tmaj",&thrust_maj_max,"tmaj/D");
  thrust_tree->Branch("tmin",&thrust_min_max,"tmin/D");
  */

  //define instream
  //string input_file = "Text_Files/med5_jewel.txt";
  string input_file = "Text_Files/pythia_files_120k.txt";
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
  //for(int ifile = 10; ifile < 11; ifile++){
  for(int ifile = startfile; ifile < endfile; ifile++){

    string str = filename[ifile];
    file = TFile::Open(str.c_str());

    if (debug) cout << "\n **** =========================== New File ================================= **** \n ";
    //cout << "File Address: " << str << endl;
    cout << "File Name: " << filename[ifile] << endl;
    cout << "File Number: " << ifile << "/" << startfile << " to " << endfile << endl;
    //cout << (Form("dijet/t%d"),radius) << endl; 

    //define trees and file ===========================
    //t = (TTree*)file->Get(Form("dijet/t%d"),radius);
    t = (TTree*)file->Get("dijet/t3");
    nt = (TTree*)file->Get("dijet/nt");
    hi = (TTree*)file->Get("ana/hi");

    //Set branches of the tree 
    t->SetBranchAddress("jtpt", &pt);
    t->SetBranchAddress("jteta", &eta);
    t->SetBranchAddress("jtphi", &phi);
    t->SetBranchAddress("nref", &nref);
    hi->SetBranchAddress("vz", &vz);
    nt->SetBranchAddress("pthat", &pThat);

    t->AddFriend(hi);
    t->AddFriend(nt);
    
    Long64_t nentries = t->GetEntries();
    //nentries = 10000;
 
    cout << "Events in File: " << nentries << endl;
    eventCount = 0;

    //Calculate weights for this file ==================
    //weightQ = xsQCD[ifile]/nentries;
    weightSTAT = xsSTAT[ifile]/nentries;
    
    //event loop
    for(Long64_t nentry = 0; nentry<nentries; ++nentry){

      if(nentry%10000 == 0) cout << nentry << endl;
      
      t->GetEntry(nentry);
      hi->GetEntry(nentry);
      nt->GetEntry(nentry);
      
      jetCount = 0;
      isGoodEvt = 0; 

      // get the pt vector for each event which passes you jet selections based on pT and eta
      vector <float> pt_v;
      vector <float> eta_v;
      vector <float> phi_v;

      for(int ij = 0; ij<nref; ++ij){
	if(pt[ij] > ptCutLead && fabs(eta[ij]) < etaCut){
	  isGoodEvt = 1; 
	}
	if(pt[ij] > ptCut && fabs(eta[ij]) < etaCut){
	  pt_v.push_back(pt[ij]);	
	  eta_v.push_back(eta[ij]);	
	  phi_v.push_back(phi[ij]);
	}
      }

      if (debug) cout<<"Total Number of Jets    : "<<nref<<endl;
      if (debug) cout<<"Number of Selected Jets : "<<pt_v.size()<<endl;

      int NJets_Sel = pt_v.size();

      if(NJets_Sel < 2) {
	if (debug) cout<<"This event had only 1 Jet"<<endl;
	continue;	
      }
      if(!isGoodEvt){continue;}

      if(debug) cout<< " \n ******* New Event ******** " << endl;
      if(debug) cout<< " ******* " << nentry << " ******** " << endl;
      
      //reset maximum values
      eventCount++;
      evt++;
      thrust_max = 0;

      vector <double> px;
      vector <double> py;
      vector <double> pz;
      
      
      for(Long64_t naxis = 0; naxis < NJets_Sel; ++naxis){
	float axis_jet_pt = pt_v[naxis];
	float axis_jet_eta = eta_v[naxis];
	float axis_jet_phi = phi_v[naxis];
      	px.push_back((double)axis_jet_pt * TMath::Cos(axis_jet_phi));
	py.push_back((double)axis_jet_pt * TMath::Sin(axis_jet_phi));
	pz.push_back((double)axis_jet_pt * TMath::SinH(axis_jet_eta));
      }

      //PART 1 ====================================================
      //Runs through all the jets in an event, checking them to see if they are the axis that maximizes thrust
      //max axis finding loop
      for(Long64_t naxis = 0; naxis < NJets_Sel; ++naxis){

	float axis_jet_pt = pt_v[naxis];
	float axis_jet_eta = eta_v[naxis];
	float axis_jet_phi = phi_v[naxis];
	
	h_pT->Fill(axis_jet_pt,weightSTAT);
	
	if(debug) cout<< " \n --------- New Test Axis (Thrust)--------- " << endl; 
	//reset values for this particular event
	thrust_temp = 0;  // maj_temp = 0;   min_temp = 0;
	
	//calculates axis for this particular jet
	TVector3 nT (px[naxis], py[naxis], pz[naxis]);
	if(debug) cout<<"Test Axis Unnormed = {" << nT(0) << ", " << nT(1) << ", " << nT(2)<< "}" << endl;
	nT = Norm(nT);
	
	if(debug) cout<<"Test Axis = {" << nT(0) << ", " << nT(1) << ", " << nT(2)<< "}" << endl;
	
	//resets for next jet loop
	dot = 0;   mag = 0;

	//PART 2 ====================================================
	//Loops through all the jets to find the thrust value for the chosen axis 
	//jet loop
	for(Long64_t njet = 0; njet < NJets_Sel; ++njet){
	  
	  float jet_pt = pt_v[njet];
	  float jet_eta = eta_v[njet];
	  float jet_phi = phi_v[njet];
	
	  if(debug) cout<< " \n --------- New Jet (Thrust)--------- " << endl; 
	  if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<jet_pt<<"\n \t eta = "<<jet_eta<<"\n \t phi = "<<jet_phi<<endl;
	
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

	  //fill jet variables
	  h_eta->Fill(eta_v[njet], weightSTAT);
	  h_phi->Fill(phi_v[njet], weightSTAT);
	  h_TBadpT->Fill(pt_v[njet], weightSTAT);
	  
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
	  max_nref = naxis; 
	}
	if(debug) cout<< "max thrust = " << thrust_max << endl;
	
      }//end axis loop
      
      if (debug) cout << "FINAL THRUST VALUE: " << thrust_max << endl; 
      
      // FILL BAD THRUST VALUES TO DEBUG =============================
      if(thrust_max < 0.5) {

	cout << "Thrust: " << thrust_max << " Event Number: " << nentry << endl; 
	isBadThrust = 1;
	badCount++;
	
	//h_TBadpT->Fill(pt_[max_nref], weightSTAT);
	h_TBad->Fill(thrust_max, weightSTAT);
	//eta and phi of bad jet axes
	//h_etaBad->Fill(eta_v[max_nref], weightSTAT);
	//h_phiBad->Fill(phi_v[max_nref], weightSTAT);
	h_weightBad->Fill(weightSTAT);
	h_pthatBad->Fill(pThat);
	h_nrefBad->Fill(nref);
	h_jetCountBad->Fill(NJets_Sel);
	h_fileNum->Fill(ifile);
	
	if (debug) cout << "______________________________" << endl; 
	if (debug) cout << "| XXX : THRUST LESS THAN 0.5 |" << endl;
	if (debug) cout << "|  Max Thrust: " << thrust_max << endl;
	//cout << "|  Max Thrust: " << thrust_max << endl;
	if (debug) cout << "______________________________" << endl; 
      }
      
      //PART 3 ====================================================
      //Begin code to select the Thrust Major and Minor axes
      
      //define the plane perpendicular to this axis in order to calculate Tmaj and Tmin
      Plane* perp = new Plane(max_thrust_axis);
      
      //reset maximum values for new axis test
      thrust_maj_max = 0;   thrust_min_max = 0; 
      
      //Thrust maj axis loop
      for(Long64_t naxis = 0; naxis < NJets_Sel; ++naxis){

	// if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)){ continue;}
	if(debug) cout<< " \n --------- New Test Axis (Min/Maj)--------- " << endl; 
       
	//define momentum three vector
	TVector3 p3 (px[naxis], py[naxis], pz[naxis]);
	if(debug) cout<<"Jet Axis UnNormed = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	p3 = Norm(p3);
	if(debug) cout<<"Jet Axis Normed = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	
	//define maj_axis and min_axis 
	TVector3 maj_axis = perp->Projection((p3));
	if(debug) cout<<"Maj Axis UnNormed = {" << maj_axis(0) << ", " << maj_axis(1) << ", " << maj_axis(2)<< "}" << endl;
	maj_axis = Norm(maj_axis);
	TVector3 min_axis = max_thrust_axis.Cross(maj_axis);
	min_axis = Norm(min_axis); 
	
	if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	if(debug) cout<<"Maj Axis = {" << maj_axis(0) << ", " << maj_axis(1) << ", " << maj_axis(2)<< "}" << endl;
	if(debug) cout<<"Min Axis = {" << min_axis(0) << ", " << min_axis(1) << ", " << min_axis(2)<< "}\n" << endl;
	
	//reset for new axis test
	dot_maj = 0;   dot_min = 0;   mag = 0;

	//PART 4 ====================================================
	//Test the axis defined by the above loop to determine if this axis is the maximum
	//jet loop
	for(Long64_t njet = 0; njet < NJets_Sel; ++njet){
	  
	  //make a ptcut
	  float jet_pt = pt_v[njet];
	  float jet_eta = eta_v[njet];
	  float jet_phi = phi_v[njet];
	  
	  if(debug) cout<< " \n --------- New Jet (Maj/Min)--------- " << endl; 
	  if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<jet_pt<<"\n \t eta = "<<jet_eta<<"\n \t phi = "<<jet_phi<<endl;
	  
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
	  /*
	  if(maj_temp > 0.5) {
	    h_TmajBadpT->Fill(pt[naxis]);
	    h_TmajBad->Fill(maj_temp);
	  } 
	  if(min_temp > 0.5) {
	    h_TminBadpT->Fill(pt[naxis]);
	    h_TminBad->Fill(maj_temp);
	  } 
	  */
	}   
      }//end of major/minor axis loop
      
      //fill all the maximum values before finishing
      h_thrust->Fill(thrust_max);
      h_min->Fill(thrust_min_max);
      h_maj->Fill(thrust_maj_max);
      h_nref->Fill(nref);
      h_jetCount->Fill(NJets_Sel);
      h_weight->Fill(pThat,weightSTAT);
      h_unweight->Fill(pThat);
	
      if(debug) {
	if (thrust_max < 0.5)     {  cout << "FLAG_thrust1: " << thrust_max <<  " , " << jetCount << endl; }
	if (thrust_maj_max > 0.5) {  cout << "FLAG_maj: " << thrust_maj_max <<  " , " << jetCount << endl; }
	if (thrust_min_max > 0.5) {  cout << "FLAG_min: " << thrust_min_max <<  " , " << jetCount << endl; }
      }	

      pt_v.clear();
      eta_v.clear();
      phi_v.clear();

      px.clear();
      py.clear();
      pz.clear();

      //thrust_tree->Fill(); 
      
    }//end of event loop
    
    gROOT->GetListOfFiles()->Remove(file);
    
    cout << "Events Selected: " << eventCount << endl;
    cout << "File Finished" << endl; 
    
  }//end file loop

  //scale the histograms
  Double_t integral;
  
  // integral = h_thrust->Integral();
  // h_thrust->Scale(1/integral);
  h_thrust->Scale(1./h_thrust->Integral());
  h_maj->Scale(1./h_maj->Integral());
  h_min->Scale(1./h_min->Integral());

  // integral = h_maj->Integral(); 
  // h_maj->Scale(1/integral);

  // integral = h_min->Integral();
  // h_min->Scale(1/integral);

  //integral = h_80->Integral();
  //h_80->Scale(integral);

  //integral = h_60->Integral();
  //h_60->Scale(integral);

  // integral = h_40->Integral();
  //h_40->Scale(integral); 
 
  //Create the plot for Thrust vs. dN/dT
  //define histograms
  TH1F * h_T = DivideByBinWidth(h_thrust, "thrust_scaled");
  TH1F * h_Tmaj = DivideByBinWidth(h_maj, "thrust_maj_scaled");
  TH1F * h_Tmin = DivideByBinWidth(h_min, "thrust_min_scaled");
  
  //Float_t divEntries = 1./(h_thrust->GetBinWidth(1));
  
  //h_T->Scale(divEntries);
  //h_Tmaj->Scale(divEntries);
  //h_Tmin->Scale(divEntries);
  
  //h_40 = DivideByBinWidth(h_40, "thrust_40_new");   h_40->Scale(divEntries);
  //h_60 = DivideByBinWidth(h_60, "thrust_60_new");   h_60->Scale(divEntries);
  //h_80 = DivideByBinWidth(h_80, "thrust_80_new");   h_80->Scale(divEntries);

  save_File->cd();

  //thrust_tree->Write();

  h_T->Print("base");
  h_Tmaj->Print("base");
  h_Tmin->Print("base");
  h_pT->Print("base");
  h_pTcut->Print("base");
  h_nref->Print("base");
  h_jetCount->Print("base");
  h_eta->Print("base");
  h_phi->Print("base");
  h_weight->Print("base"); 
  
  h_T->Write();
  h_Tmaj->Write();
  h_Tmin->Write();
  h_thrust->Write();
  h_maj->Write();
  h_min->Write();
  h_pT->Write();
  h_pTcut->Write();
  h_nref->Write();
  h_jetCount->Write();
  h_eta->Write();
  h_phi->Write();
  h_weight->Write();
  h_unweight->Write();
  
  save_File->Write();
  save_File->Close();

  cout<<endl<<endl<<endl<<endl<<"GOING TO WRITE BAD OUTPUT FILE"<<endl<<endl<<endl<<endl;
  
  //if(debug) {

    bad_File->cd();
    
    h_TBad->Print("base");
    h_TmajBad->Print("base");
    h_TminBad->Print("base");
    h_TBadpT->Print("base");
    h_TmajBadpT->Print("base");
    h_TminBadpT->Print("base");
    h_etaBad->Print("base");
    h_phiBad->Print("base");
    h_nrefBad->Print("base");
    h_jetCountBad->Print("base");
    h_pthatBad->Print("base");
    h_weightBad->Print("base");
    h_fileNum->Print("base"); 
    
    h_TBad->Write();
    h_TmajBad->Write();
    h_TminBad->Write();
    h_TBadpT->Write();
    h_TmajBadpT->Write();
    h_TminBadpT->Write();
    h_etaBad->Write();
    h_phiBad->Write();
    h_nrefBad->Write();
    h_jetCountBad->Write();
    h_pthatBad->Write();
    h_weightBad->Write();
    h_fileNum->Write();
    
    bad_File->Write();
    bad_File->Close();
    // }

  timer.Stop();
  cout<<"Bad Percent: "<<badCount<<"/"<<evt<<endl;
  cout<<"Macro finished: "<<endl;
  
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
    
}//end of plot thrust
