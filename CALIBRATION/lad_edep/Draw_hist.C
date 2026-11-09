#include <TCanvas.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLine.h>
#include <TProfile.h>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

const char spec_prefix = 'P'; // Spectrometer to replay

const int N_PLANES                 = 5;
const int N_PADDLES                = 11;
const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const int N_SIDES                  = 2;
const string side_names[N_SIDES]   = {"Top", "Btm"};

const double fit_func_xmin = 2.2; // Minimum x value for fit function
const double fit_func_xmax = 6.0; // Maximum x value for fit function

const double quad_fit_param[2] = {5.31637, 12.9101};
double cut_line_pt1[2]   = {2, 120}; // Point 1 for cut line
const double cut_line_pt2[2]   = {6, 0};   // Point 2 for cut line

const double true_adc_MeV = 80; // True ADC value in MeV

void Draw_hist(){
	TString inputFileName = "lad_edep_plots_new_timing_H.root";
	TFile *file = TFile::Open(inputFileName, "READ");
	if (!file || file->IsZombie()) {
		cout << "Error: Cannot open ROOT file." << endl;
		return;
	}
	TFile *outfile = TFile::Open("hist.root", "RECREATE");
	if (!outfile || outfile->IsZombie()) {
		printf("Error: Cannot create output ROOT file.\n");
		file->Close();
		return;
	}
	string hist_names[2] = {"KIN/TDC_DIFF_VS_ADC_AMP/c_TDC_DIFF_VS_ADC_AMP_plane_%s",
		"KIN/TDC_DIFF_VS_ADC_INT/c_TDC_DIFF_VS_ADC_INT_plane_%s"};
	double raw_adc_max[2][N_PLANES][N_PADDLES] = {0};
	for (int i_hist = 0; i_hist < 1; ++i_hist) {
		string hist_name = hist_names[i_hist];
		for (int i_plane = 0; i_plane < N_PLANES; ++i_plane) {

			TCanvas *canvas                     = (TCanvas *)file->Get(Form(hist_name.c_str(), plane_names[i_plane].c_str()));
			TList *primitives                   = canvas->GetListOfPrimitives();
			TH2D *hist_arr[N_PLANES][N_PADDLES] = {nullptr};

			TH2D *hist = nullptr;
			int count  = 0;
			for (TObject *obj : *primitives) {
				if (obj->InheritsFrom("TPad")) {
					TPad *pad            = (TPad *)obj;
					TList *padPrimitives = pad->GetListOfPrimitives();
					for (TObject *padObj : *padPrimitives) {
						if (padObj->InheritsFrom("TH2")) {
							hist_arr[0][count++] = (TH2D *)padObj;
							// Now you can use 'hist' as before
						}
					}
				}
			}
      cout<< "Plane: " << plane_names[i_plane] << ", Histograms found: " << count << endl;
      TCanvas* c1 = new TCanvas(Form("c_%s_plane_%s", hist_name.c_str(), plane_names[i_plane].c_str()),
                              Form("%s Plane %s", hist_name.c_str(), plane_names[i_plane].c_str()), 800, 600);
      c1->Divide(4, 3);
	  for (int i_paddle=0; i_paddle < N_PADDLES; ++i_paddle) {
		c1->cd(i_paddle+1);
		if (hist_arr[0][i_paddle]) {
	  	TGraph *graph = new TGraph();
        hist_arr[0][i_paddle]->DrawClone("COLZ");
		hist_arr[0][i_paddle]->Smooth();
		int nPoints=0;
		for (int i= 1; i <= hist_arr[0][i_paddle]->GetNbinsX(); ++i) {
			double x = hist_arr[0][i_paddle]->GetXaxis()->GetBinCenter(i);
			TH1D* proj = hist_arr[0][i_paddle]->ProjectionY(Form("h_%s_plane_%s_bar_%i_bin_%i", hist_name.c_str(), plane_names[i_plane].c_str(),i_paddle,i), i, i, "e");
        	//proj->Smooth();
			TSpectrum *spec = new TSpectrum();
			spec->Search(proj, 1, "goff", 0.1);
			if (spec->GetNPeaks() >= 1) {
				//Get the highest peak
				double peak = spec->GetPositionX()[0];
				double peak_height = spec->GetPositionY()[0];
				for (int j = 1; j < spec->GetNPeaks(); ++j) {
					if (spec->GetPositionY()[j] > peak_height) {
						peak = spec->GetPositionX()[j];
						peak_height = spec->GetPositionY()[j];
					}
				}
				graph->SetPoint(nPoints++, x, peak);
				cout << "Peak found at: " << peak << endl;
			}else{
				cout << "No peak found for bin " << i << endl;
			}
			delete proj;
			delete spec;
		}
		if(nPoints>0){
			graph->SetName(Form("graph_%s_plane_%s_bar_%i", hist_name.c_str(), plane_names[i_plane].c_str(), i_paddle));
			graph->SetLineColor(kRed);
			graph->SetMarkerStyle(20);
			graph->SetMarkerColor(kRed);
			graph->Draw("P");
			}
		}
	  }
	  c1->Update();
			
      outfile->cd();
      c1->Write();
   } 
  }

}