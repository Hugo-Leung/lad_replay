void plot_gem_trk_rdataframe(){
    const int MAX_DATA = 1500; // Adjust as needed 
    TString inList = "filelist_23367_23370.txt";
    TString outName = "gem_trk_23367_23370.root";
    TString spec[] = {"H","P"}; // Spectrometer name (H or P)
    const int nSpecs = 2;
    std::ifstream fileList(inList);
    if (!fileList.is_open()) {
        std::cerr << "Error: Unable to open file list: " << inList << std::endl;
        return;
    }
    TChain* T = new TChain("T");
    std::string line;
    while (std::getline(fileList, line)) {
        if (!line.empty()) {
            T->Add(line.c_str());
        }
    }
    fileList.close();
    TFile *outputFile = new TFile(outName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
        return;
    }
    TTreeReader reader(T);
    // Define variables
    TTreeReaderValue<Double_t>* triggerBits_float = new TTreeReaderValue<Double_t>(reader, "g.trigbits");
    TTreeReaderValue<Int_t>* nTrack[2];
    TTreeReaderValue<Int_t>* nGoodHits[2];
    TTreeReaderArray<Double_t>* trk_d0[2];
    TTreeReaderArray<Double_t>* trk_dt[2];
    TTreeReaderArray<Double_t>* trk_t[2];
    TTreeReaderArray<Double_t>* trk_chisq[2];
    TTreeReaderArray<Double_t>* trk_theta[2];
    TTreeReaderArray<Double_t>* trk_phi[2];
    TTreeReaderArray<Double_t>* trk_projz[2];
    TTreeReaderArray<Double_t>* trk_projx[2];
    TTreeReaderArray<Double_t>* trk_projy[2];
    TTreeReaderArray<Double_t>* trk_x1_local[2];
    TTreeReaderArray<Double_t>* trk_y1_local[2];
    TTreeReaderArray<Double_t>* trk_x2_local[2];
    TTreeReaderArray<Double_t>* trk_y2_local[2];
    TTreeReaderArray<Double_t>* trk_has_hodo_hit[2];
    TTreeReaderArray<Double_t>* trk_is_good[2];
    TTreeReaderValue<Double_t>* react_x[2];
    TTreeReaderValue<Double_t>* react_y[2];
    TTreeReaderValue<Double_t>* react_z[2];
    TTreeReaderArray<Double_t>* hitTime[2][2];
    TTreeReaderArray<Double_t>* hitTrkID[2][2];
    TTreeReaderArray<Double_t>* trk_adc[2][2];
    for (int i=0; i<2; i++){
        nTrack[i] = new TTreeReaderValue<Int_t>(reader, Form("Ndata.%s.gem.trk.d0", spec[i].Data()));
        nGoodHits[i] = new TTreeReaderValue<Int_t>(reader, Form("Ndata.%s.ladhod.goodhit_trackid", spec[i].Data()));
        trk_d0[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.d0", spec[i].Data()));
        trk_dt[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.dt", spec[i].Data()));
        trk_t[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.t", spec[i].Data()));
        trk_chisq[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.chisq", spec[i].Data()));
        trk_theta[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.theta", spec[i].Data()));
        trk_phi[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.phi", spec[i].Data()));
        trk_projz[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.projz", spec[i].Data()));
        trk_projx[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.projx", spec[i].Data()));
        trk_projy[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.projy", spec[i].Data()));
        trk_x1_local[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.x1_local", spec[i].Data()));
        trk_y1_local[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.y1_local", spec[i].Data()));
        trk_x2_local[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.x2_local", spec[i].Data()));
        trk_y2_local[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.y2_local", spec[i].Data()));
        trk_has_hodo_hit[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.has_hodo_hit", spec[i].Data()));
        trk_is_good[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.is_good", spec[i].Data()));
        react_x[i] = new TTreeReaderValue<Double_t>(reader, Form("%s.react.x", spec[i].Data()));
        react_y[i] = new TTreeReaderValue<Double_t>(reader, Form("%s.react.y", spec[i].Data()));
        react_z[i] = new TTreeReaderValue<Double_t>(reader, Form("%s.react.z", spec[i].Data()));
        hitTime[i][0] = new TTreeReaderArray<Double_t>(reader, Form("%s.ladhod.goodhit_hittime_0", spec[i].Data()));
        hitTime[i][1] = new TTreeReaderArray<Double_t>(reader, Form("%s.ladhod.goodhit_hittime_1", spec[i].Data()));
        hitTrkID[i][0] = new TTreeReaderArray<Double_t>(reader, Form("%s.ladhod.goodhit_trackid", spec[i].Data()));
        trk_adc[i][0] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.adc1", spec[i].Data()));
        trk_adc[i][1] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.adc2", spec[i].Data()));
    }
    //define histograms
    const int nHisto = 4; //all good tracks, no hodo hit, 1 hodo hit, 2 hodo hits
    TH1D* h_trk_d0[nHisto];
    TH1D* h_trk_dt[nHisto];
    TH1D* h_trk_t[nHisto];
    TH1D* h_trk_chisq[nHisto];
    TH1D* h_trk_theta[nHisto];
    TH1D* h_trk_phi[nHisto];
    TH1D* h_trk_projz[nHisto];
    TH1D* h_trk_is_good[nHisto];
    TH2D* h_vertex_z[nHisto];
    TH2D* h_z_phi[nHisto];
    TH1D* h_hit_time[nHisto][2]; // hit time for hodo hits in plane 0 and plane 1
    TH2D* h_hit_local[nHisto][2]; // local hit position for gem hits in near and far
    TH1D* h_hit_adc[nHisto][2]; // adc for gem hits in near and far
    TH2D* h_hit_adc_ylocal[nHisto][2]; // adc for gem hits in near and far
    TH1D* h_trk_has_hodo_hit = new TH1D("h_trk_has_hodo_hit", "Track Has Hodo Hit;Has Hodo Hit;Counts", 3, -0.5, 2.5);
    TH1D* h_bad_chisq = new TH1D ("h_bad_chisq","Failed tracking chisq", 6, -6.5,-0.5);
    for(int i=0; i<nHisto; i++){
        h_trk_d0[i] = new TH1D(Form("h_trk_d0_%d", i), Form("Track D0 (Category %d);D0 (cm);Counts", i), 100, 0, 50);
        h_trk_dt[i] = new TH1D(Form("h_trk_dt_%d", i), Form("Track Dt (Category %d);Dt (ns);Counts", i), 100, -100, 100);
        h_trk_t[i] = new TH1D(Form("h_trk_t_%d", i), Form("Track T (Category %d);T (ns);Counts", i), 100, 0, 200);
        h_trk_chisq[i] = new TH1D(Form("h_trk_chisq_%d", i), Form("Track Chi2 (Category %d);Chi2;Counts", i), 100, 0, 30);
        h_trk_theta[i] = new TH1D(Form("h_trk_theta_%d", i), Form("Track Theta (Category %d);Theta (degrees);Counts", i), 100, 70, 170);
        h_trk_phi[i] = new TH1D(Form("h_trk_phi_%d", i), Form("Track Phi (Category %d);Phi (degrees);Counts", i), 100, -35, 35);
        h_trk_projz[i] = new TH1D(Form("h_trk_projz_%d", i), Form("Track ProjZ (Category %d);ProjZ (cm);Counts", i), 100, -25, 25);
        h_trk_is_good[i] = new TH1D(Form("h_trk_is_good_%d", i), Form("Track Is Good (Category %d);Is Good;Counts", i), 3, -0.5, 2.5);
        h_vertex_z[i] = new TH2D(Form("h_vertex_z_%d", i), Form("Vertex Z (Category %d);Track ProjZ (cm);Reaction Z (cm)", i), 100, -25, 25, 100, -25, 25);
        h_z_phi[i] = new TH2D(Form("h_z_phi_%d", i), Form("Z vs Phi (Category %d);Z (cm);Phi (degrees)", i), 30, -15, 15, 70, -35, 35);
        for (int j=0; j<2; j++){
            h_hit_time[i][j] = new TH1D(Form("h_hit_time_%d_%d", i, j), Form("Hit Time (Category %d, Plane %d);Hit Time (ns);Counts", i, j), 100, 1500, 2000);
            h_hit_local[i][j] = new TH2D(Form("h_hit_local_%d_%d", i, j), Form("GEM Local Hit Position (Category %d, Plane %d);X (cm);Y (cm)", i, j), 1000, -70, 70, 1000, -40, 40);
            h_hit_adc[i][j] = new TH1D(Form("h_hit_adc_%d_%d", i, j), Form("GEM ADC (Category %d, Plane %d);ADC;Counts", i, j), 100, 0, 10000);
            h_hit_adc_ylocal[i][j] = new TH2D(Form("h_hit_adc_ylocal_%d_%d", i, j), Form("GEM ADC vs Y Local (Category %d, Plane %d);Y Local (cm);ADC", i, j), 40, -40, 40, 100, 0, 10000);
        }
    }
    std::cout << "Starting event loop..." << std::endl;

    int j =0;

    while (reader.Next()) {
        if (j++ % 1000 == 0) {
            std::cout << "Processed " << j << " entries.\r" << std::flush;
        }
        //if (j>100000) break; // Limit to 500,000 entries for testing
        UInt_t triggerBits = static_cast<UInt_t>(**triggerBits_float);
        bool hmsEreal = (triggerBits & 8) != 0;
        bool shmsEreal = (triggerBits & 2) != 0;
        int iSpec = -1;
        if (hmsEreal && !shmsEreal) {
            iSpec = 0; // HMS
        } else if (!hmsEreal && shmsEreal) {
            iSpec = 1; // SHMS
        } else {
            std::cerr << "Warning: Event does not match expected trigger conditions. Skipping event." << std::endl;
            std::cout << "Trigger bits: " << (triggerBits & 0xff) << std::endl;
            continue; // Skip events that don't match the trigger conditions
        }
        int nTracks_iSpec = **nTrack[iSpec];
        for (int i = 0; i < nTracks_iSpec; i++) {
            int iCategory = 0; // Default category for all good tracks
            if ((*trk_has_hodo_hit[iSpec])[i] == 0) {
                iCategory = 1; // No hodo hit
            } else if ((*trk_has_hodo_hit[iSpec])[i] == 1) {
                iCategory = 2; // 1 hodo hit
            } else if ((*trk_has_hodo_hit[iSpec])[i] >= 2) {
                iCategory = 3; // 2 or more hodo hits
            }
            Double_t hitTime_tmp[2] = {-1, -1}; // Temporary array to store hit times for the current track
            for (int j=0; j<**nGoodHits[iSpec]; j++){
                for (int plane=0; plane<2; plane++){
                    if ((*hitTrkID[iSpec][0])[j]==i){
                        hitTime_tmp[plane]=(*hitTime[iSpec][plane])[j];
                    }
                }
            }
            if((*trk_chisq[iSpec])[i]<0){
                h_bad_chisq->Fill((*trk_chisq[iSpec])[i]);
            }

            h_trk_is_good[0]->Fill((*trk_is_good[iSpec])[i]);
            if ((*trk_is_good[iSpec])[i] != 1) {
                continue; // Skip tracks that are not marked as good
            }
            h_trk_d0[0]->Fill((*trk_d0[iSpec])[i]);
            h_trk_d0[iCategory]->Fill((*trk_d0[iSpec])[i]);
            h_trk_dt[0]->Fill((*trk_dt[iSpec])[i]);
            h_trk_dt[iCategory]->Fill((*trk_dt[iSpec])[i]);
            h_trk_t[0]->Fill((*trk_t[iSpec])[i]);
            h_trk_t[iCategory]->Fill((*trk_t[iSpec])[i]);
            h_trk_chisq[0]->Fill((*trk_chisq[iSpec])[i]);
            h_trk_chisq[iCategory]->Fill((*trk_chisq[iSpec])[i]);
            h_trk_theta[0]->Fill((*trk_theta[iSpec])[i]*TMath::RadToDeg());
            h_trk_phi[0]->Fill((*trk_phi[iSpec])[i]*TMath::RadToDeg());
            h_trk_projz[0]->Fill((*trk_projz[iSpec])[i]);
            h_hit_local[0][0]->Fill((*trk_x1_local[iSpec])[i], (*trk_y1_local[iSpec])[i]);
            h_hit_local[0][1]->Fill((*trk_x2_local[iSpec])[i], (*trk_y2_local[iSpec])[i]);
            h_z_phi[0]->Fill((*trk_projz[iSpec])[i], (*trk_phi[iSpec])[i]*TMath::RadToDeg());
            h_hit_adc[0][0]->Fill((*trk_adc[iSpec][0])[i]);
            h_hit_adc[0][1]->Fill((*trk_adc[iSpec][1])[i]);
            h_hit_adc_ylocal[0][0]->Fill((*trk_y1_local[iSpec])[i], (*trk_adc[iSpec][0])[i]);
            h_hit_adc_ylocal[0][1]->Fill((*trk_y2_local[iSpec])[i], (*trk_adc[iSpec][1])[i]);
            //if(true){
            if(((iCategory==1&&(*trk_chisq[iSpec])[i]<1) || (iCategory==2&&(*trk_chisq[iSpec])[i]<5) || (iCategory==3&&(*trk_chisq[iSpec])[i]<10))&&((*trk_projz[iSpec])[i]>-19 && (*trk_projz[iSpec])[i]<19)){
                h_trk_theta[iCategory]->Fill((*trk_theta[iSpec])[i]*TMath::RadToDeg());
                h_trk_phi[iCategory]->Fill((*trk_phi[iSpec])[i]*TMath::RadToDeg());
                h_trk_projz[iCategory]->Fill((*trk_projz[iSpec])[i]);
                h_vertex_z[iCategory]->Fill((*trk_projz[iSpec])[i], (**react_z[iSpec]));
                h_hit_local[iCategory][0]->Fill((*trk_x1_local[iSpec])[i], (*trk_y1_local[iSpec])[i]);
                h_hit_local[iCategory][1]->Fill((*trk_x2_local[iSpec])[i], (*trk_y2_local[iSpec])[i]);
                h_z_phi[iCategory]->Fill((*trk_projz[iSpec])[i], (*trk_phi[iSpec])[i]*TMath::RadToDeg());
                if (hitTime_tmp[0]>0&& hitTime_tmp[0]<1e6){
                    h_hit_time[0][0]->Fill(hitTime_tmp[0]);
                    h_hit_time[iCategory][0]->Fill(hitTime_tmp[0]);
                }
                if (hitTime_tmp[1]>0&& hitTime_tmp[1]<1e6){
                    h_hit_time[0][1]->Fill(hitTime_tmp[1]);
                    h_hit_time[iCategory][1]->Fill(hitTime_tmp[1]);
                }
                h_hit_adc[iCategory][0]->Fill((*trk_adc[iSpec][0])[i]);
                h_hit_adc[iCategory][1]->Fill((*trk_adc[iSpec][1])[i]);
                h_hit_adc_ylocal[iCategory][0]->Fill((*trk_y1_local[iSpec])[i], (*trk_adc[iSpec][0])[i]);
                h_hit_adc_ylocal[iCategory][1]->Fill((*trk_y2_local[iSpec])[i], (*trk_adc[iSpec][1])[i]);
            }
            h_trk_has_hodo_hit->Fill((*trk_has_hodo_hit[iSpec])[i]);
            h_vertex_z[0]->Fill((*trk_projz[iSpec])[i], (**react_z[iSpec]));
        }
        
        
    }
    cout<<"\nFinished processing entries."<<std::endl;
    cout<<"Total good tracks: "<<h_trk_chisq[0]->GetEntries()<<std::endl;

    double fixedz[] = {-10, 0, 10};
    double z_range[] = {-15, -5, 5, 15}; // Define the z ranges for the fixed positions, should be one more than the number of fixed z positions
    const int nFixedz = 0;

    /*double fixedz[]={-5,5};//optics2;
    double z_range[] = {-10,0,10};
    const int nFixedz =2; */
    for (int i=0; i<nHisto; i++){
        cout<<"Fitting track projection Z for category "<<i<<std::endl;
        double means[nFixedz];
        double means_err[nFixedz];
        double sigmas[nFixedz];
        double sigmas_err[nFixedz];
        TPaveText *pt = new TPaveText(.5, .55, .95, .65, "NB NDC");
        pt->SetFillColor(0);
        pt->SetBorderSize(0);
        pt->SetFillStyle(0);
        for (int ifoil=0; ifoil<nFixedz; ifoil++){
            TFitResultPtr fitResult = h_trk_projz[i]->Fit("gaus", "QLISNR", "", fixedz[ifoil]-2.5, fixedz[ifoil]+2.5);
            double mean = fitResult->Parameter(1);
            double sigma = fitResult->Parameter(2);
            cout<<"Z = "<<fixedz[ifoil]<<": Mean = "<<mean<<" cm, Sigma = "<<sigma<<" cm"<<std::endl;
            double fitRangeMin = TMath::Max (mean - 2*sigma,-25.0);
            double fitRangeMax = TMath::Min (mean + 2*sigma, 25.0); 
            TFitResultPtr fitResult2 = h_trk_projz[i]->Fit("gaus", "QLIS+R", "", fitRangeMin, fitRangeMax);
            means[ifoil] = fitResult2->Parameter(1);
            means_err[ifoil] = fitResult2->ParError(1);
            sigmas[ifoil] = fitResult2->Parameter(2);
            sigmas_err[ifoil] = fitResult2->ParError(2);
            pt->AddText(Form("Z = %.1f cm: Mean = %.2f #pm %.2f cm, Sigma = %.2f #pm %.2f cm", fixedz[ifoil], means[ifoil], means_err[ifoil], sigmas[ifoil], sigmas_err[ifoil]));
            TLine *line = new TLine(fixedz[ifoil], 0, fixedz[ifoil], h_trk_projz[i]->GetMaximum());
            line->SetLineColor(kGray+1);
            h_trk_projz[i]->GetListOfFunctions()->Add(line);
            TLine *line2 = new TLine(-25, fixedz[ifoil], 25, fixedz[ifoil]);
            line2->SetLineColor(kGray+1);
            h_vertex_z[i]->GetListOfFunctions()->Add(line2);
            TLine *line3 = new TLine(fixedz[ifoil], -25, fixedz[ifoil], 25);
            line3->SetLineColor(kGray+1);
            h_vertex_z[i]->GetListOfFunctions()->Add(line3);
        }
        h_trk_projz[i]->GetListOfFunctions()->Add(pt);
        
        double totalCounts =0; 
        double signalCounts=0;
        TBox* backgroundBox = new TBox(z_range[0], z_range[0], z_range[nFixedz], z_range[nFixedz]);
        backgroundBox->SetLineColor(kGreen+2);
        backgroundBox->SetLineWidth(2);
        backgroundBox->SetFillStyle(0);
        h_vertex_z[i]->GetListOfFunctions()->Add(backgroundBox);
        for (int iFoil=0; iFoil<nFixedz; iFoil++){
            for (int jFoil=0; jFoil<nFixedz; jFoil++){
                int binX1 = h_vertex_z[i]->GetXaxis()->FindBin(z_range[iFoil]);
                int binX2 = h_vertex_z[i]->GetXaxis()->FindBin(z_range[iFoil+1]);
                int binY1 = h_vertex_z[i]->GetYaxis()->FindBin(z_range[jFoil]);
                int binY2 = h_vertex_z[i]->GetYaxis()->FindBin(z_range[jFoil+1]);
                double content = h_vertex_z[i]->Integral(binX1, binX2, binY1, binY2);
                totalCounts+=content;
                if (iFoil==jFoil){
                    signalCounts+=content;
                }
                std::cout << "Signal counts in region (" << z_range[iFoil] << ", " << z_range[iFoil+1] << ") x (" << z_range[jFoil] << ", " << z_range[jFoil+1] << "): " << content << std::endl;
                
            }
            TBox* box = new TBox(z_range[iFoil],z_range[iFoil], z_range[iFoil+1], z_range[iFoil+1]);
            box->SetLineColor(kRed);
            box->SetLineWidth(2);
            box->SetFillStyle(0);
            h_vertex_z[i]->GetListOfFunctions()->Add(box);
        }
        std::cout << "Total counts in vertex distribution: " << totalCounts << std::endl;
        std::cout << "Estimated signal counts in defined regions: " << signalCounts << std::endl;
        std::cout << "Estimated background counts: " << totalCounts - signalCounts << std::endl;
        double ratio=signalCounts / (totalCounts - signalCounts) ;
        std::cout << "Estimated signal-to-background ratio: " << ratio << std::endl;
        TPaveText *pt_ratio = new TPaveText(.4, .85, .75, .90, "NB NDC");
        pt_ratio->SetFillColor(0);
        pt_ratio->SetBorderSize(0);
        pt_ratio->SetFillStyle(0);
        pt_ratio->AddText(Form("Signal-to-background ratio: %.2f", ratio));
        h_vertex_z[i]->GetListOfFunctions()->Add(pt_ratio);

        for (int j=0; j<2; j++){
            double h_hit_time_mean = h_hit_time[i][j]->GetMean();
            TLine *line_mean = new TLine(h_hit_time_mean, 0, h_hit_time_mean, h_hit_time[i][j]->GetMaximum());
            TLine *line_lo = new TLine(h_hit_time_mean-24*3, 0, h_hit_time_mean-24*3, h_hit_time[i][j]->GetMaximum());
            TLine *line_hi = new TLine(h_hit_time_mean+24*3, 0, h_hit_time_mean+24*3, h_hit_time[i][j]->GetMaximum());
            h_hit_time[i][j]->GetListOfFunctions()->Add(line_mean);
            h_hit_time[i][j]->GetListOfFunctions()->Add(line_lo);
            h_hit_time[i][j]->GetListOfFunctions()->Add(line_hi);
        }
    }


    outputFile->Write();
    outputFile->Close();
    delete triggerBits_float;
    for (int i=0; i<2; i++){
        delete nTrack[i];
        delete nGoodHits[i];
        delete trk_d0[i];
        delete trk_dt[i];
        delete trk_t[i];
        delete trk_chisq[i];
        delete trk_theta[i];
        delete trk_phi[i];
        delete trk_projz[i];
        delete trk_projx[i];
        delete trk_projy[i];
        delete trk_x1_local[i];
        delete trk_y1_local[i];
        delete trk_x2_local[i];
        delete trk_y2_local[i];
        delete trk_has_hodo_hit[i];
        delete trk_is_good[i];
        delete react_x[i];
        delete react_y[i];
        delete react_z[i];
        delete hitTrkID[i][0];
        for (int j=0; j<2; j++){
            delete hitTime[i][j];
            delete trk_adc[i][j];
        }
    }
    delete T;
}