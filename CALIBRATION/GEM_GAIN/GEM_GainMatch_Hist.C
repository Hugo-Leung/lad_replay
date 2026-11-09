vector<double> AsymALL, dAsymALL;
vector<int> APVX_asymALL, APVY_asymALL;
vector<vector<double> > weightX, weightY;
vector<int> countX,countY;
//Number of possible APVs in X and Y directions. Looks like in beam transport coordinates
//These variables get overridden by the "strip config" file:


const int nAPVmaxX = 12;
const int nAPVmaxY = 24;

double width_cut =6./7; // only fit around the peak, needs to find a better way to define this value.

void GEM_GainMatch_Hist(const TString inName="GEM_gainmatch_temp.root", const TString outName="GEM_gainMatch_histfit.root"){
	TF1* asymFitFunc = new TF1("asymFitFunc","gaus(0)+pol1(3)", -1.0, 1.0);
	
    const int nmodules=2; //number of GEM modules per spectrometer
    
    TFile *inputFile = new TFile(inName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open the input ROOT file!" << std::endl;
        return;
    }
    TFile *outputFile = new TFile(outName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
        return;
    }

    for (int iModule=0; iModule<nmodules; iModule++){
        //Fit the module-average asymmetry histograms first
        TString hname_mod;
        hname_mod.Form("BeforeGainMatch/h_ADCasym_clust_mod%d", iModule);
        TH1D* h_ADCasym_clust = (TH1D*)inputFile->Get(hname_mod.Data());
        if (!h_ADCasym_clust) {
            std::cerr << "Error: Cannot find histogram " << hname_mod.Data() << " in the input file!" << std::endl;
            continue;
        }

        if(h_ADCasym_clust->GetEntries() < 100) {
            std::cerr << "Warning: Histogram " << hname_mod.Data() << " has less than 100 entries. Skipping fit." << std::endl;
            outputFile->cd();
            h_ADCasym_clust->Write(); // Still write the histogram to the output file
            continue;
        }
        int imaxbin = h_ADCasym_clust->GetMaximumBin();
        double maxAsym = h_ADCasym_clust->GetXaxis()->GetBinCenter(imaxbin);
        double maxCounts = h_ADCasym_clust->GetBinContent(imaxbin);
        int ibinlo=imaxbin;
        int ibinhi=imaxbin;
        while (ibinhi<=h_ADCasym_clust->GetNbinsX() && h_ADCasym_clust->GetBinContent(ibinhi)/maxCounts > width_cut){
            ibinhi++;
        }
        while (ibinlo>=1 && h_ADCasym_clust->GetBinContent(ibinlo)/maxCounts > width_cut){
            ibinlo--;
        }
        TFitResultPtr fitresult = h_ADCasym_clust->Fit("gaus","SLR","",h_ADCasym_clust->GetXaxis()->GetBinCenter(ibinlo), h_ADCasym_clust->GetXaxis()->GetBinCenter(ibinhi));
        asymFitFunc->SetParameters(fitresult->Parameter(0), fitresult->Parameter(1), fitresult->Parameter(2), 0, 0);
        double ylo= h_ADCasym_clust->GetBinContent( h_ADCasym_clust->FindBin(-0.65) );
        double yhi= h_ADCasym_clust->GetBinContent( h_ADCasym_clust->FindBin(0.65) );
        double slope = (yhi - ylo)/(0.65 - (-0.65));
        double intercept = ylo - slope*(-0.65);
        asymFitFunc->SetParameter(3, intercept);
        asymFitFunc->SetParameter(4, slope);
        TFitResultPtr fitresult2 = h_ADCasym_clust->Fit(asymFitFunc,"SLR","", -0.65, 0.65);
        outputFile->cd();
        h_ADCasym_clust->Write();
        double fittedMean = fitresult2->Parameter(1);
        double fittedSigma = fitresult2->Parameter(2);

        for (int iAPVx=0; iAPVx<nAPVmaxX; iAPVx++){
            for (int iAPVy=0; iAPVy<nAPVmaxY; iAPVy++){
                TString hname_xy;
                hname_xy.Form("DebugHistograms/h_ADCasym_APVXY_mod%d_apvx%d_apvy%d", iModule, iAPVx, iAPVy);
                TH1D* h_ADCasym_APVXY = (TH1D*)inputFile->Get(hname_xy.Data());
                if (!h_ADCasym_APVXY) {
                    std::cerr << "Error: Cannot find histogram " << hname_xy.Data() << " in the input file!" << std::endl;
                    continue;
                }
                if(h_ADCasym_APVXY->GetEntries() < 100) {
                    std::cerr << "Warning: Histogram " << hname_xy.Data() << " has less than 100 entries. Skipping fit." << std::endl;
                    outputFile->cd();
                    h_ADCasym_APVXY->Write(); // Still write the histogram to the output file
                    continue;
                }
                

		        int imaxbin = h_ADCasym_APVXY->GetMaximumBin();
		        double maxAsym = h_ADCasym_APVXY->GetXaxis()->GetBinCenter(imaxbin);
		        double maxCounts = h_ADCasym_APVXY->GetBinContent(imaxbin);
		        int ibinlo=imaxbin;
		        int ibinhi=imaxbin;
		        while (ibinhi<=h_ADCasym_APVXY->GetNbinsX() && h_ADCasym_APVXY->GetBinContent(ibinhi)/maxCounts > width_cut){
			        ibinhi++;
		        }
		        while (ibinlo>=1 && h_ADCasym_APVXY->GetBinContent(ibinlo)/maxCounts > width_cut){
			        ibinlo--;
		        }
                TFitResultPtr fitresult = h_ADCasym_APVXY->Fit("gaus","SLR","",h_ADCasym_APVXY->GetXaxis()->GetBinCenter(ibinlo), h_ADCasym_APVXY->GetXaxis()->GetBinCenter(ibinhi));
                double MPV = fitresult->Parameter(1);
                double Sigma = fitresult->Parameter(2);
                if(MPV < -0.5 || MPV > 0.1 || TMath::Abs(Sigma) < 0.01 || TMath::Abs(Sigma) > 0.2){
                    MPV=fittedMean;
                    Sigma=fittedSigma;
                }
		        asymFitFunc->SetParameters(fitresult->Parameter(0), MPV, Sigma, 0, 0);
		        double ylo= h_ADCasym_APVXY->GetBinContent( h_ADCasym_APVXY->FindBin(-0.67) );
		        double yhi= h_ADCasym_APVXY->GetBinContent( h_ADCasym_APVXY->FindBin(0.67) );
		        double slope = (yhi - ylo)/(0.67 - (-0.67));
		        double intercept = ylo - slope*(-0.67);
		        asymFitFunc->SetParameter(3, intercept);
		        asymFitFunc->SetParameter(4, slope);
		        TFitResultPtr fitresult2 = h_ADCasym_APVXY->Fit(asymFitFunc,"SLR","", -0.67, 0.67);

                //write the histogram and fit result to the output file
                outputFile->cd();
                h_ADCasym_APVXY->Write();
            }
        }
    }
    inputFile->Close();
    outputFile->Write();
    outputFile->Close();

}