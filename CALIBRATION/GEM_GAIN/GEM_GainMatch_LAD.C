#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TSystem.h"
#include <iostream>
#include <fstream> 
#include "TROOT.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVectorD.h"
#include "TMatrixD.h"
//#include "Math/Functor.h"
//#include "TMinuit.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCut.h"
#include "TTreeFormula.h"


//Global variables used by the chi2 function:
vector<double> AsymALL, dAsymALL;
vector<int> APVX_asymALL, APVY_asymALL;
vector<vector<double> > weightX, weightY;
vector<int> countX,countY;
//Number of possible APVs in X and Y directions. Looks like in beam transport coordinates
//These variables get overridden by the "strip config" file:
const int nAPVmaxX = 12;
const int nAPVmaxY = 24;


void chi2_FCN( int &npar, double *gin, double &f, double *par, int flag ){

	//cout << "calculating chi2..." << endl;
	//why do we assume this for the start parameter?
	double chi2 = 0.0;
	for( int i=0; i<AsymALL.size(); i++ ){
		//Define loop variables iy and ix as y and x APV asymmetry respectively?
		// this is the combination of X and Y APV cards for the i'th entry in the asymmetry array
		int iy = APVY_asymALL[i];
		int ix = APVX_asymALL[i];
		//What is this variable for? this variable becomes true if certain quality checks are passed and this particular combination should be included in the fit
		bool useasym = false;
		// Gain starting values? why do we pick 1.0 for start? 1 is default
		double Gx=1.0, Gy=1.0;
		//what is the purpose of each conditional statement?
		if( iy >= 0 && iy < nAPVmaxY && ix >= 0 && ix < nAPVmaxX ){
			//Need to understand par before can understand this variable assignment
			Gx = par[ix + nAPVmaxY];
			Gy = par[iy];
			//why does this variable change? because if we are here then we have a good asymmetry for this combination of X and Y APVs
			useasym = true;
		} else if( iy >= 0 && iy < nAPVmaxY ){ //this is a module-average asymmetry for Y: using weighted average over all X APVs for Gx:
			double sum_Gx = 0.0;
			double sum_Wx = 0.0;
			for( int j = 0; j<nAPVmaxX; j++ ){
				sum_Gx += weightX[iy][j] * par[j+nAPVmaxY];
				sum_Wx += weightX[iy][j];
			}
			//Which definition is actually be used. What is the mechanic here?

			//This is replacing Gx with the statistics-weighted average of all the X APVs:
			Gx = sum_Gx/sum_Wx;
			Gy = par[iy];

			//only use the module-average asymmetry in the fit if there were no X APVs with at least 100 events on this Y APV:
			// 
			// if( countX[iy] == 0 )

			//Let's see what happens if we use the module-average asymmetry here:
			if( countX[iy] == 0 ) useasym = true;

		} else if( ix >= 0 && ix < nAPVmaxX ){ //this is a module-average asymmetry for X: using weighted average over all Y APVs for Gy:
			double sum_Gy = 0.0;
			double sum_Wy = 0.0;
			for( int j=0; j<nAPVmaxY; j++ ){
				sum_Gy += weightY[ix][j] * par[j];
				sum_Wy += weightY[ix][j];
			}

			Gx = par[ix+nAPVmaxY];
			Gy = sum_Gy/sum_Wy;

			if( countY[ix] == 0 ) useasym = true;

		}
		//The actual asymmetry  
		double Atheory = (Gx - Gy)/(Gx + Gy);
		//it is the sum over all measured asymmetries of the difference between the measured and "expected" asymmetries, which are funcitons of the APV gain coeffs:
		if( useasym ) chi2 += pow( (AsymALL[i] - Atheory)/dAsymALL[i], 2 );

	}

	//cout << "chi2 = " << chi2 << endl;

	f = chi2;
}

void GEM_GainMatch_LAD( const char *configfilename, const char *outfname="GEM_gainmatch_temp.root" , int debug=1 ){
	ifstream configfile(configfilename);

	if( !configfile ) return;

	TChain *C = new TChain("T");

	TString currentline;

	while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
		if( !currentline.BeginsWith("#") ){
			C->Add(currentline);
		}
	}
	configfile.close();

	TFile *outputFile = new TFile(outfname, "RECREATE");
	if (!outputFile || outputFile->IsZombie()) {
		std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
		return; 
	}

	UInt_t MAXNHITS=100000;
	UInt_t MAXNTRACKS=10000;

	double ADCcut = 1000.0;
	double ccor_cut = 0.95;
	double deltat_cut = 2.5;


	TString spec[] = {"H","P"}; // Spectrometer name (H or P)
	UInt_t triggerMask[]={8,2};// HMS Ereal, SHMS Ereal
	const int nSpecs = 2;
	const int nmodules=2; //number of GEM modules per spectrometer
	Double_t triggerBits_float=0;

	//create a bunch of vector<double>, presumably to store info?
	std::vector<Double_t> sp_clustID1[2];
	std::vector<Double_t> sp_clustID2[2];
	Int_t nSp[2];

	std::vector<double> clust_module[2];
	std::vector<double> clust_layer[2];
	std::vector<double> clust_axis[2];
	std::vector<Double_t> clust_ID[2];
	std::vector<Double_t> clust_striplo[2];
	std::vector<Double_t> clust_striphi[2];
	std::vector<Double_t> clust_nstrip[2];
	std::vector<Double_t> clust_stripmax[2];
	Int_t nClust[2];

	Double_t module_nClust[2][2]; //spec, layer
	Double_t module_n2Dhits[2][2]; //spec, layer

	std::vector<Double_t> adc[2];
	std::vector<Double_t> gain[2];
	std::vector<Double_t> sp_adcAsym[2];
	std::vector<Double_t> sp_adcMean[2];
	std::vector<Double_t> sp_dt[2];
	std::vector<Double_t> sp_ccor[2];
	std::vector<Double_t> sp_id[2];

	std::vector<double> sp_layer[2];
	std::vector<Double_t> sp_goodhit[2];

	Int_t ntracks[2];
	std::vector<Double_t> trk_spid_0[2];
	std::vector<Double_t> trk_spid_1[2];
	std::vector<Double_t> trk_sp_0u[2];
	std::vector<Double_t> trk_sp_0v[2];
	std::vector<Double_t> trk_sp_1u[2];
	std::vector<Double_t> trk_sp_1v[2];

	Double_t hit_nlayer[2];
	std::vector<TString> varnames;

	TH1D* h_triggerBits;
	TH1D* h_ClustID[2][2];
	TH1D* h_ClustAxis[2][2];
	TH1D* h_nGoodSp[2];
	TH1D* h_ADCasym_clust[2];//module-wise histograms of ADC asymmetry from clusters
	TH1D* h_ADCasym_sp[2];//module-wise histograms of ADC asymmetry from space points
	TH1D* h_ADCasym_diff[2];//module-wise histograms of difference in ADC asymmetry between clusters and space points
	TH1D* h_ADCmean_clust[2];//module-wise histograms of ADC mean from clusters
	TH1D* h_ADCmean_all=new TH1D("h_ADCmean_all","ADC Mean from all modules;ADC Mean;Counts",1500,0,30000);

	// Bunch of debugging historgrams
	outputFile->mkdir("BeforeGainMatch");
	outputFile->cd("BeforeGainMatch");

	for(int i=0; i<nmodules;i++){
		h_ADCasym_clust[i] = new TH1D(Form("h_ADCasym_clust_mod%d", i), Form("ADC Asymmetry from clusters Module %d;ADC Asymmetry;Counts", i), 250, -1.01, 1.01);
		h_ADCasym_sp[i] = new TH1D(Form("h_ADCasym_sp_mod%d", i), Form("ADC Asymmetry from space points Module %d;ADC Asymmetry;Counts", i), 250, -1.01, 1.01);
		h_ADCasym_diff[i] = new TH1D(Form("h_ADCasym_diff_mod%d", i), Form("Difference in ADC Asymmetry Module %d;ADC Asymmetry Difference;Counts", i), 250, -1.01, 1.01);
		h_ADCmean_clust[i] = new TH1D(Form("h_ADCmean_clust_mod%d", i), Form("ADC Mean from clusters Module %d;ADC Mean;Counts", i), 1500, 0, 30000);
		h_ADCasym_clust[i]->Sumw2();
		h_ADCasym_sp[i]->Sumw2();
		h_ADCmean_clust[i]->Sumw2();
	}

	outputFile->mkdir("AfterGainMatch");
	if (debug){
		outputFile->cd();
		outputFile->mkdir("DebugHistograms");
		h_triggerBits = new TH1D("h_triggerBits", "Trigger Bits;Trigger Bits;Counts", 3, -0.5, 2.5);
		for(int i =0; i<nSpecs; i++){
			h_ClustID[i][0] = new TH1D(Form("h_ClustID1_%s", spec[i].Data()), Form("Cluster ID1 %s;Cluster ID1;Counts", spec[i].Data()), 1000, -0.5, 1000.5);
			h_ClustAxis[i][0] = new TH1D(Form("h_ClustAxis1_%s", spec[i].Data()), Form("Cluster Axis1 %s;Cluster Axis1;Counts", spec[i].Data()), 2, -0.5, 1.5);
			h_ClustID[i][1] = new TH1D(Form("h_ClustID2_%s", spec[i].Data()), Form("Cluster ID2 %s;Cluster ID2;Counts", spec[i].Data()), 1000, -0.5, 1000.5);
			h_ClustAxis[i][1] = new TH1D(Form("h_ClustAxis2_%s", spec[i].Data()), Form("Cluster Axis2 %s;Cluster Axis2;Counts", spec[i].Data()), 2, -0.5, 1.5);
			h_nGoodSp[i] = new TH1D(Form("h_nGoodSp_%s", spec[i].Data()), Form("Number of good space points %s;N good space points;Counts", spec[i].Data()), 2000, -0.5, 2000.5);
		}
	}



	C->SetBranchAddress("g.trigbits", &triggerBits_float);
	varnames.push_back("g.trigbits");  
	for (int iSpec=0; iSpec<2; iSpec++){
		nSp[iSpec]=-1;
		sp_clustID1[iSpec]=std::vector<Double_t>(MAXNHITS);
		sp_clustID2[iSpec]=std::vector<Double_t>(MAXNHITS);
		sp_layer[iSpec]=std::vector<double>(MAXNHITS);
		sp_dt[iSpec]=std::vector<Double_t>(MAXNHITS);
		sp_ccor[iSpec]=std::vector<Double_t>(MAXNHITS);
		sp_id[iSpec]=std::vector<Double_t>(MAXNHITS);

		clust_module[iSpec]=std::vector<double>(MAXNHITS);
		clust_layer[iSpec]=std::vector<double>(MAXNHITS);
		clust_axis[iSpec]=std::vector<double>(MAXNHITS);
		clust_ID[iSpec]=std::vector<Double_t>(MAXNHITS); 
		clust_striplo[iSpec]=std::vector<Double_t>(MAXNHITS);
		clust_striphi[iSpec]=std::vector<Double_t>(MAXNHITS);
		clust_nstrip[iSpec]=std::vector<Double_t>(MAXNHITS);
		clust_stripmax[iSpec]=std::vector<Double_t>(MAXNHITS);
		nClust[iSpec]=-1;


		sp_clustID1[iSpec].resize(MAXNHITS);
		sp_clustID2[iSpec].resize(MAXNHITS);
		clust_module[iSpec].resize(MAXNHITS);
		clust_layer[iSpec].resize(MAXNHITS);
		clust_axis[iSpec].resize(MAXNHITS);
		sp_layer[iSpec].resize(MAXNHITS);

		adc[iSpec]=std::vector<Double_t>(MAXNHITS);
		gain[iSpec]=std::vector<Double_t>(MAXNHITS);
		sp_adcAsym[iSpec]=std::vector<Double_t>(MAXNHITS);
		sp_adcMean[iSpec]=std::vector<Double_t>(MAXNHITS);
		sp_goodhit[iSpec]=std::vector<Double_t>(MAXNHITS);

		ntracks[iSpec]=-1;
		trk_spid_0[iSpec]=std::vector<Double_t>(MAXNHITS);
		trk_spid_1[iSpec]=std::vector<Double_t>(MAXNHITS);
		trk_sp_0u[iSpec]=std::vector<Double_t>(MAXNHITS);
		trk_sp_0v[iSpec]=std::vector<Double_t>(MAXNHITS);
		trk_sp_1u[iSpec]=std::vector<Double_t>(MAXNHITS);
		trk_sp_1v[iSpec]=std::vector<Double_t>(MAXNHITS);

		varnames.push_back(TString::Format("Ndata.%s.gem.sp.layer",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.clusID1",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.clusID2",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.layer",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.spID", spec[iSpec].Data()));

		varnames.push_back(TString::Format("Ndata.%s.gem.clust.axis",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.module",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.layer",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.axis",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.index",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.striplo",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.striphi",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.nstrip",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.maxstrip",spec[iSpec].Data()));

		varnames.push_back(TString::Format("%s.gem.clust.adc",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.clust.apvGain",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.asym",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.adc",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.dt",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.corrcoeff",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.sp.isgoodhit",spec[iSpec].Data()));

		varnames.push_back(TString::Format("%s.gem.hit.nlayer",spec[iSpec].Data()));

		varnames.push_back(TString::Format("Ndata.%s.gem.trk.id",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.trk.spID_0",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.trk.spID_1",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.trk.spID_0u",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.trk.spID_0v",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.trk.spID_1u",spec[iSpec].Data()));
		varnames.push_back(TString::Format("%s.gem.trk.spID_1v",spec[iSpec].Data()));


		for (int i=0; i<nmodules; i++){
			varnames.push_back(TString::Format("%s.gem.m%d.nclusters",spec[iSpec].Data(),i));
			varnames.push_back(TString::Format("%s.gem.m%d.n2Dhits",spec[iSpec].Data(),i));
		}
	}

	C->SetBranchStatus("*",0); //disable all branches first
	//enable only the branches we need:
	for( int i=0; i<varnames.size(); i++ ){
		C->SetBranchStatus( varnames[i], 1 );
	}
	cout << "Setting branch addresses: ";
	for (int iSpec=0; iSpec<nSpecs; iSpec++){
		C->SetBranchAddress( TString::Format("Ndata.%s.gem.sp.layer",spec[iSpec].Data()), &nSp[iSpec] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.layer",spec[iSpec].Data()), &sp_layer[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.clusID1",spec[iSpec].Data()), &sp_clustID1[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.clusID2",spec[iSpec].Data()), &sp_clustID2[iSpec][0] );

		C->SetBranchAddress( TString::Format("Ndata.%s.gem.clust.axis",spec[iSpec].Data()), &nClust[iSpec] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.module",spec[iSpec].Data()), &clust_module[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.layer",spec[iSpec].Data()), &clust_layer[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.axis",spec[iSpec].Data()), &clust_axis[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.index",spec[iSpec].Data()), &clust_ID[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.striplo",spec[iSpec].Data()), &clust_striplo[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.striphi",spec[iSpec].Data()), &clust_striphi[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.nstrip",spec[iSpec].Data()), &clust_nstrip[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.maxstrip",spec[iSpec].Data()), &clust_stripmax[iSpec][0] );

		C->SetBranchAddress( TString::Format("%s.gem.clust.adc",spec[iSpec].Data()), &adc[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.clust.apvGain",spec[iSpec].Data()), &gain[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.asym",spec[iSpec].Data()), &sp_adcAsym[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.adc",spec[iSpec].Data()), &sp_adcMean[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.dt",spec[iSpec].Data()), &sp_dt[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.corrcoeff",spec[iSpec].Data()), &sp_ccor[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.isgoodhit",spec[iSpec].Data()), &sp_goodhit[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.sp.spID",spec[iSpec].Data()), &sp_id[iSpec][0]);

		C->SetBranchAddress( TString::Format("%s.gem.hit.nlayer",spec[iSpec].Data()), &hit_nlayer[iSpec] );

		C->SetBranchAddress( TString::Format("Ndata.%s.gem.trk.id",spec[iSpec].Data()), &ntracks[iSpec] );
		C->SetBranchAddress( TString::Format("%s.gem.trk.spID_0",spec[iSpec].Data()), &trk_spid_0[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.trk.spID_1",spec[iSpec].Data()), &trk_spid_1[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.trk.spID_0u",spec[iSpec].Data()), &trk_sp_0u[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.trk.spID_0v",spec[iSpec].Data()), &trk_sp_0v[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.trk.spID_1u",spec[iSpec].Data()), &trk_sp_1u[iSpec][0] );
		C->SetBranchAddress( TString::Format("%s.gem.trk.spID_1v",spec[iSpec].Data()), &trk_sp_1v[iSpec][0] );

		for (int i=0; i<nmodules; i++){
			C->SetBranchAddress( TString::Format("%s.gem.m%d.nclusters",spec[iSpec].Data(),i), &module_nClust[iSpec][i] );
			C->SetBranchAddress( TString::Format("%s.gem.m%d.n2Dhits",spec[iSpec].Data(),i), &module_n2Dhits[iSpec][i] );
			
		}
	}
	cout << "done." << endl;

	//Define array to hold asymetry values for each apv

	TH1D* h_ADCasym_APVXY[nAPVmaxX*nAPVmaxY*nmodules];
	TH1D* h_ADCasym_APVX[nAPVmaxX*nmodules];
	TH1D* h_ADCasym_APVY[nAPVmaxY*nmodules];

	for (int imodule=0; imodule<nmodules; imodule++){
		for (int ix=0; ix<nAPVmaxX; ix++){
			for (int iy=0; iy<nAPVmaxY; iy++){
				TString hname_xy;
				hname_xy.Form("h_ADCasym_APVXY_mod%d_apvx%d_apvy%d", imodule, ix, iy);
				h_ADCasym_APVXY[iy + nAPVmaxY*ix + nAPVmaxX*nAPVmaxY*imodule] = new TH1D(hname_xy.Data(), hname_xy.Data(), 250, -1.01, 1.01);
				h_ADCasym_APVXY[iy + nAPVmaxY*ix + nAPVmaxX*nAPVmaxY*imodule]->Sumw2();
			}
		}
		for (int ix=0; ix<nAPVmaxX; ix++){
			TString hname_x;
			hname_x.Form("h_ADCasym_APVX_mod%d_apvx%d", imodule, ix);
			h_ADCasym_APVX[ix + nAPVmaxX*imodule] = new TH1D(hname_x.Data(), hname_x.Data(), 250, -1.01, 1.01);
			h_ADCasym_APVX[ix + nAPVmaxX*imodule]->Sumw2();
		}
		for (int iy=0; iy<nAPVmaxY; iy++){
			TString hname_y;
			hname_y.Form("h_ADCasym_APVY_mod%d_apvy%d", imodule, iy);
			h_ADCasym_APVY[iy + nAPVmaxY*imodule] = new TH1D(hname_y.Data(), hname_y.Data(), 250, -1.01, 1.01);
			h_ADCasym_APVY[iy + nAPVmaxY*imodule]->Sumw2();
		}
	}



	// int nEvents = C->GetEntries();
	//cout << "Total number of events: " << nEvents << endl;
	TRandom3 randgen(1); //initialize random number generator with seed 1
	double debugFraction = 0.05; //process only 5% of events in debug mode
	//if (debug) cout << "Debug mode: processing only " << debugFraction*100 << "\% of events." << endl;
	int iEvent=0;
	TF1* asymFitFunc = new TF1("asymFitFunc","gaus(0)+pol1(3)", -1.0, 1.0);
	while(C->GetEntry(iEvent++)){
		if ( iEvent % 10000 == 0 ) cout << "\rProcessing " << iEvent << " events." << std::flush;
		//if (debug && randgen.Uniform() > debugFraction) continue; //skip event in debug mode based on random selection
		//if (debug && iEvent>5e4) break; //limit number of events in debug mode
		UInt_t trigbits = static_cast<UInt_t>(triggerBits_float);
		int iSpec = -1;
		if( (trigbits & triggerMask[0]) != 0 ){ //HMS
			iSpec = 0;
		} else if( (trigbits & triggerMask[1]) != 0 ){ //SHMS
			iSpec = 1;
		}
		if(debug){
			h_triggerBits->Fill(iSpec);
			if(iSpec>=0) h_nGoodSp[iSpec]->Fill(nSp[iSpec]);
		}
		if( iSpec == -1 ){ 
			cout<<"No valid spectrometer triggered at event "<<iEvent<<", trigbits = "<<trigbits<<endl;
			continue; //no valid spectrometer triggered
		}
		if (nSp[iSpec]<0 || nSp[iSpec]>MAXNHITS){
			cout << "Invalid number of space points at event " << iEvent << ": nSp = " << nSp[iSpec] << endl;
			continue;
		}

		for (int iSp=0; iSp<nSp[iSpec]; iSp++){
			if (iSp<0||iSp>nSp[iSpec]){
				cout<<"invalid spID: "<<iSp<<"; nSp: "<<nSp[iSpec]<<endl;
				continue;
			}
			int clustID1 = (int) sp_clustID1[iSpec][iSp];
			int clustID2 = (int) sp_clustID2[iSpec][iSp];
			int module_index = (int) sp_layer[iSpec][iSp];
			
			if(module_index==0&&(sp_id[iSpec][iSp]>nSp[iSpec] || sp_id[iSpec][iSp]!=iSp)){
				cout<<"module: "<<module_index<<" Sp.id: "<<sp_id[iSpec][iSp]<<"; iSp: "<<iSp<<"; nSp: "<<nSp[iSpec]<<endl;
			}
			if(module_index==1&&(sp_id[iSpec][iSp]>(nSp[iSpec]-module_n2Dhits[iSpec][0]) || sp_id[iSpec][iSp]!=(iSp-module_n2Dhits[iSpec][0]))){
				cout<<"module: "<<module_index<<" Sp.id: "<<sp_id[iSpec][iSp]<<"; iSp: "<<iSp<<"; nSp: "<<nSp[iSpec]<<endl;
			}


			if(sp_layer[iSpec][iSp]>0){
				clustID1+=module_nClust[iSpec][0];
				clustID2+=module_nClust[iSpec][0];
			}

			double ADCasym_sp = sp_adcAsym[iSpec][iSp];;
			double ADCmean_sp = sp_adcMean[iSpec][iSp];;
			double ADC_clust1 = adc[iSpec][clustID1];;
			double ADC_clust2 = adc[iSpec][clustID2];;
			double gain_clust1 = gain[iSpec][clustID1];;
			double gain_clust2 = gain[iSpec][clustID2];;

			double ADCasym_clust = (ADC_clust1 - ADC_clust2) / (ADC_clust1 + ADC_clust2);
			double ADCasym_clust_noGain = (ADC_clust1/gain_clust1 - ADC_clust2/gain_clust2) / (ADC_clust1/gain_clust1 + ADC_clust2/gain_clust2);
			double ADCasym_diff = ADCasym_clust - ADCasym_sp;
			double ADCmean_clust = 0.5 * (ADC_clust1 + ADC_clust2);
			double ADCmean_clust_noGain = 0.5 * (ADC_clust1/gain_clust1 + ADC_clust2/gain_clust2);

			int iAPVx_lo = (int)clust_striplo[iSpec][clustID1]/128;
			int iAPVx_hi = (int)clust_striphi[iSpec][clustID1]/128;
			int iAPVx_max = (int)clust_stripmax[iSpec][clustID1]/128;

			int iAPVy_lo = (int)clust_striplo[iSpec][clustID2]/128;
			int iAPVy_hi = (int)clust_striphi[iSpec][clustID2]/128;
			int iAPVy_max = (int)clust_stripmax[iSpec][clustID2]/128;

			//ADC and nstrip cuts can be added here if needed
			if (ADCmean_clust<ADCcut || ADC_clust1<ADCcut || ADC_clust2<ADCcut){
				continue;
			}

			if (clust_nstrip[iSpec][clustID1]<2 || clust_nstrip[iSpec][clustID2]<2){
				continue;
			}
			if (sp_ccor[iSpec][iSp]<ccor_cut){
				continue;
			}
			if (TMath::Abs(sp_dt[iSpec][iSp])>deltat_cut){
				continue;
			}

			//
			if(TMath::Abs(hit_nlayer[iSpec]-2)>1e-3|| sp_goodhit[iSpec][iSp]<0.5){
				continue;
			}

			//only use hits where all stips are within same APV card:
			if (iAPVx_lo != iAPVx_hi || iAPVy_lo != iAPVy_hi){
				continue;
			}
			if (iAPVx_lo != iAPVx_max || iAPVy_lo != iAPVy_max){
				//this should never happen
				std::cout<< " Warning: cluster max strip not within low and high strips at event " << iEvent << ", index " << iSp << std::endl;
				continue;
			}

			if (iAPVx_lo<0 || iAPVx_lo>=nAPVmaxX || iAPVy_lo<0 || iAPVy_lo>=nAPVmaxY){
				cout<< " Invalid APV index at event " << iEvent << ", index " << iSp << ": iAPVx_lo= " << iAPVx_lo << ", iAPVy_lo= " << iAPVy_lo << endl;
				continue;
			}

			h_ADCasym_clust[module_index]->Fill(ADCasym_clust);
			h_ADCasym_sp[module_index]->Fill(ADCasym_sp);
			h_ADCasym_diff[module_index]->Fill(ADCasym_diff);
			h_ADCmean_clust[module_index]->Fill(ADCmean_clust);
			h_ADCasym_APVXY[iAPVy_lo + nAPVmaxY*iAPVx_lo + nAPVmaxX*nAPVmaxY*module_index]->Fill(ADCasym_clust_noGain);
			h_ADCasym_APVX[iAPVx_lo + nAPVmaxX*module_index]->Fill(ADCasym_clust_noGain);
			h_ADCasym_APVY[iAPVy_lo + nAPVmaxY*module_index]->Fill(ADCasym_clust_noGain); 
			if(debug){
				h_ClustID[iSpec][0]->Fill(clust_ID[iSpec][clustID1]);
				h_ClustAxis[iSpec][0]->Fill(clust_axis[iSpec][clustID1]);
				h_ClustID[iSpec][1]->Fill(clust_ID[iSpec][clustID2]);
				h_ClustAxis[iSpec][1]->Fill(clust_axis[iSpec][clustID2]);
			}
			h_ADCmean_all->Fill(ADCmean_clust_noGain);
		}
	}
	TFitResultPtr fitresult = h_ADCmean_all->Fit("landau","QSR","",2000,25000);
	double mpv = fitresult->Parameter(1);





	if(debug){
		outputFile->cd();
		outputFile->cd("DebugHistograms");
		h_triggerBits->Write();
		for(int i =0; i<nSpecs; i++){
			h_ClustID[i][0]->Write();
			h_ClustAxis[i][0]->Write();
			h_ClustID[i][1]->Write();
			h_ClustAxis[i][1]->Write();
			h_nGoodSp[i]->Write();
		}
		for (int iAPVx=0; iAPVx<nAPVmaxX; iAPVx++){
			for (int iAPVy=0; iAPVy<nAPVmaxY; iAPVy++){
				for(int imodule=0; imodule<nmodules; imodule++){
		  			h_ADCasym_APVXY[iAPVy + nAPVmaxY*iAPVx + nAPVmaxX*nAPVmaxY*imodule]->Write();
		  		}
		  	}
		  	for(int imodule=0; imodule<nmodules; imodule++){
		  		h_ADCasym_APVX[iAPVx + nAPVmaxX*imodule]->Write();
		  	}
		}
		for (int iAPVy=0; iAPVy<nAPVmaxY; iAPVy++){
			for(int imodule=0; imodule<nmodules; imodule++){
				h_ADCasym_APVY[iAPVy + nAPVmaxY*imodule]->Write();
			}
		}
	}

	double width_cut =6./7; // only fit around the peak, needs to find a better way to define this value.
	for (int imodule=0; imodule<nmodules; imodule++){
		int imaxbin = h_ADCasym_clust[imodule]->GetMaximumBin();
		double maxAsym = h_ADCasym_clust[imodule]->GetXaxis()->GetBinCenter(imaxbin);
		double maxCounts = h_ADCasym_clust[imodule]->GetBinContent(imaxbin);
		int ibinlo=imaxbin;
		int ibinhi=imaxbin;
		while (ibinhi<=h_ADCasym_clust[imodule]->GetNbinsX() && h_ADCasym_clust[imodule]->GetBinContent(ibinhi)/maxCounts > width_cut){
			ibinhi++;
		}
		while (ibinlo>=1 && h_ADCasym_clust[imodule]->GetBinContent(ibinlo)/maxCounts > width_cut){
			ibinlo--;
		}
		cout<<"module "<<imodule<<": ADC asymmetry peak at "<<maxAsym<<", fitting range: ["<<h_ADCasym_clust[imodule]->GetXaxis()->GetBinCenter(ibinlo)<<", "<<h_ADCasym_clust[imodule]->GetXaxis()->GetBinCenter(ibinhi)<<"]"<<endl; 
		TFitResultPtr fitresult = h_ADCasym_clust[imodule]->Fit("gaus","SLR","",h_ADCasym_clust[imodule]->GetXaxis()->GetBinCenter(ibinlo), h_ADCasym_clust[imodule]->GetXaxis()->GetBinCenter(ibinhi));
		asymFitFunc->SetParameters(fitresult->Parameter(0), fitresult->Parameter(1), fitresult->Parameter(2), 0, 0);
		double ylo= h_ADCasym_clust[imodule]->GetBinContent( h_ADCasym_clust[imodule]->FindBin(-0.65) );
		double yhi= h_ADCasym_clust[imodule]->GetBinContent( h_ADCasym_clust[imodule]->FindBin(0.65) );
		double slope = (yhi - ylo)/(0.65 - (-0.65));
		double intercept = ylo - slope*(-0.65);
		asymFitFunc->SetParameter(3, intercept);
		asymFitFunc->SetParameter(4, slope);
		TFitResultPtr fitresult2 = h_ADCasym_clust[imodule]->Fit(asymFitFunc,"SLR","", -0.65, 0.65);
		cout<<"  Final fit results: amplitude = "<<fitresult2->Parameter(0)<<", mean = "<<fitresult2->Parameter(1)<<", sigma = "<<fitresult2->Parameter(2)<<endl;
		double signal = fitresult2->Parameter(0)*fitresult2->Parameter(2)*sqrt(2*TMath::Pi());
		double background = 0.5 * fitresult2->Parameter(4) * (0.65 - (-0.65)) * (0.65 - (-0.65)) + fitresult2->Parameter(3) * (0.65 - (-0.65));
		cout<<"  Signal = "<<signal<<", Background = "<<background<<", S/B = "<<signal/background<<endl;
	}


	outputFile->cd(); 
	outputFile->cd("BeforeGainMatch");

	for(int i=0; i<nmodules;i++){
		h_ADCasym_clust[i]->Write();
		//h_ADCasym_sp[i]->Write();
		//h_ADCasym_diff[i]->Write();
		h_ADCmean_clust[i]->Write();
	}
	h_ADCmean_all->Write();
	std::cout << std::endl << "Finished processing all events." << std::endl;

}
