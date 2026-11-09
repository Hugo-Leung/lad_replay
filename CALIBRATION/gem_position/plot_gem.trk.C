void plot_gem_trk(){
    const int MAX_DATA = 1500; // Adjust as needed 
    TString inList = "filelist_23696_5sigma.txt";
    TString outName = "gem_trk_23696_5sigma.root";
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

    Double_t trk_d0[2][MAX_DATA];
    Double_t trk_d0_good[2][MAX_DATA];
    Double_t trk_projz[2][MAX_DATA], trk_projx[2][MAX_DATA], trk_projy[2][MAX_DATA];
    Double_t trk_chisq[2][MAX_DATA];
    Double_t trk_theta[2][MAX_DATA], trk_phi[2][MAX_DATA];
    Double_t trk_x1_local[2][MAX_DATA], trk_y1_local[2][MAX_DATA];
    Double_t trk_x2_local[2][MAX_DATA], trk_y2_local[2][MAX_DATA];
    Double_t trk_has_hodo_hit[2][MAX_DATA];
    Double_t trk_is_good[2][MAX_DATA];
    Double_t trk_dt[2][MAX_DATA];
    Double_t trk_t[2][MAX_DATA];
    Double_t react_x[2], react_y[2], react_z[2];
    Double_t hitTime[2][2][MAX_DATA];
    Double_t hitTrkID[2][2][MAX_DATA];
    Double_t trk_adc[2][2][MAX_DATA];
    int nTracks[2], nGoodHits[2];

    UInt_t triggerBits=0;
    Double_t triggerBits_float=0;
    UInt_t triggerMask[]={8,2};// HMS Ereal, SHMS Ereal
    T->SetBranchAddress("g.trigbits", &triggerBits_float);
    for (int iSpec=0; iSpec<nSpecs; iSpec++){
        T->SetBranchAddress(spec[iSpec]+".gem.trk.d0", &trk_d0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.d0_good", &trk_d0_good[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.projz", &trk_projz[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.projx", &trk_projx[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.projy", &trk_projy[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.chisq", &trk_chisq[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.theta", &trk_theta[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.phi", &trk_phi[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.has_hodo_hit", &trk_has_hodo_hit[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.is_good", &trk_is_good[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.dt", &trk_dt[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.t", &trk_t[iSpec]);
        T->SetBranchAddress("Ndata."+spec[iSpec]+".gem.trk.d0", &nTracks[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".react.x", &react_x[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".react.y", &react_y[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".react.z", &react_z[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.x1_local", &trk_x1_local[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.y1_local", &trk_y1_local[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.x2_local", &trk_x2_local[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.y2_local", &trk_y2_local[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.adc1", &trk_adc[iSpec][0]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.adc2", &trk_adc[iSpec][1]);

        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hittime_0", &hitTime[0][iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hittime_1", &hitTime[1][iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_trackid_0", &hitTrkID[0][iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_trackid_1", &hitTrkID[1][iSpec]);
        T->SetBranchAddress("Ndata."+spec[iSpec]+".ladkin.goodhit_hittime_0", &nGoodHits[iSpec]);
    }
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


    TH1D* h_trk_chisq_bad = new TH1D("h_trk_chisq_bad", "Track Chi2 (Bad);Chi2;Counts", 100, -10, 0);
    TH1D* h_trk_has_hodo_hit = new TH1D("h_trk_has_hodo_hit", "Track Has Hodo Hit;Has Hodo Hit;Counts", 3, -0.5, 2.5);
    TH1D* h_diff_x = new TH1D("h_diff_x","Differene between projected x and reaction x;ProjX - ReactX (cm);Counts", 100, -10.5, 10.5);
    TH1D* h_diff_y = new TH1D("h_diff_y","Differene between projected y and reaction y;ProjY - ReactY (cm);Counts", 100, -10.5, 10.5);

    int nEntries = T->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;
    for (int entry=0; entry<nEntries; entry++){
        if (entry%1000==0){
            std::cout<<"\rProcessing entry "<<entry<<"/"<<nEntries<<std::flush;
        }
        T->GetEntry(entry);
        int iSpec=-1;
        triggerBits = static_cast<UInt_t>(triggerBits_float);
        if ((triggerBits&triggerMask[0])!=0){//HMS
            iSpec=0;
        }else if ((triggerBits&triggerMask[1])!=0){//SHMS
            iSpec=1;
        }else{
            std::cout<<"no valid spectrometer triggered"<<std::endl;
            std::cout<<"event "<<entry<<": triggerBits = "<<triggerBits<<std::endl;
            continue; // no valid spectrometer triggered
        }

        for(int i=0; i<nTracks[iSpec]; i++){
            if (trk_chisq[iSpec][i]<0){
                h_trk_chisq_bad->Fill(trk_chisq[iSpec][i]);
                continue; // Skip tracks with bad chi2
            }
            int iCategory = 0; // default category for all good tracks
            int tmp_has_hodo_hit = (int)trk_has_hodo_hit[iSpec][i];
            if (tmp_has_hodo_hit>=0 && tmp_has_hodo_hit<3){
                iCategory = tmp_has_hodo_hit + 1; // categories 1, 2, 3 for hodo hit counts
            }else{
                std::cout<<"Invalid hodo hit count: "<<tmp_has_hodo_hit<<std::endl;
                continue; // Skip tracks with invalid hodo hit count
            }
            double hitTime_tmp[2]={-1, -1};
            for (int j=0; j<nGoodHits[iSpec]; j++){
                for (int plane=0; plane<2; plane++){
                    if (hitTrkID[plane][iSpec][j]==i){
                        hitTime_tmp[plane]=hitTime[plane][iSpec][j];
                    }
                }
            }

            h_trk_is_good[0]->Fill(trk_is_good[iSpec][i]);
            h_trk_is_good[iCategory]->Fill(trk_is_good[iSpec][i]);
            if (trk_is_good[iSpec][i] == 0){
                continue; // Skip tracks that are not marked as good
            }
            h_trk_d0[0]->Fill(trk_d0[iSpec][i]);
            h_trk_d0[iCategory]->Fill(trk_d0[iSpec][i]);
            h_trk_dt[0]->Fill(trk_dt[iSpec][i]);
            h_trk_dt[iCategory]->Fill(trk_dt[iSpec][i]);
            h_trk_t[0]->Fill(trk_t[iSpec][i]);
            h_trk_t[iCategory]->Fill(trk_t[iSpec][i]);
            h_trk_chisq[0]->Fill(trk_chisq[iSpec][i]);
            h_trk_chisq[iCategory]->Fill(trk_chisq[iSpec][i]);
            h_trk_theta[0]->Fill(trk_theta[iSpec][i]*TMath::RadToDeg());
            h_trk_phi[0]->Fill(trk_phi[iSpec][i]*TMath::RadToDeg());
            h_trk_projz[0]->Fill(trk_projz[iSpec][i]);
            h_hit_local[0][0]->Fill(trk_x1_local[iSpec][i], trk_y1_local[iSpec][i]);
            h_hit_local[0][1]->Fill(trk_x2_local[iSpec][i], trk_y2_local[iSpec][i]);
            h_z_phi[0]->Fill(trk_projz[iSpec][i], trk_phi[iSpec][i]*TMath::RadToDeg());
            h_hit_adc[0][0]->Fill(trk_adc[iSpec][0][i]);
            h_hit_adc[0][1]->Fill(trk_adc[iSpec][1][i]);
            h_hit_adc_ylocal[0][0]->Fill(trk_y1_local[iSpec][i], trk_adc[iSpec][0][i]);
            h_hit_adc_ylocal[0][1]->Fill(trk_y2_local[iSpec][i], trk_adc[iSpec][1][i]);
            if(((iCategory==1&&trk_chisq[iSpec][i]<1) || (iCategory==2&&trk_chisq[iSpec][i]<5) || (iCategory==3&&trk_chisq[iSpec][i]<10))&&(trk_projz[iSpec][i]>-15 && trk_projz[iSpec][i]<15)){
                h_trk_theta[iCategory]->Fill(trk_theta[iSpec][i]*TMath::RadToDeg());
                h_trk_phi[iCategory]->Fill(trk_phi[iSpec][i]*TMath::RadToDeg());
                h_trk_projz[iCategory]->Fill(trk_projz[iSpec][i]);
                h_vertex_z[iCategory]->Fill(trk_projz[iSpec][i], react_z[iSpec]);
                h_hit_local[iCategory][0]->Fill(trk_x1_local[iSpec][i], trk_y1_local[iSpec][i]);
                h_hit_local[iCategory][1]->Fill(trk_x2_local[iSpec][i], trk_y2_local[iSpec][i]);
                h_z_phi[iCategory]->Fill(trk_projz[iSpec][i], trk_phi[iSpec][i]*TMath::RadToDeg());
                if (hitTime_tmp[0]>0&& hitTime_tmp[0]<1e6){
                    h_hit_time[0][0]->Fill(hitTime_tmp[0]);
                    h_hit_time[iCategory][0]->Fill(hitTime_tmp[0]);
                }
                if (hitTime_tmp[1]>0&& hitTime_tmp[1]<1e6){
                    h_hit_time[0][1]->Fill(hitTime_tmp[1]);
                    h_hit_time[iCategory][1]->Fill(hitTime_tmp[1]);
                }
                h_hit_adc[iCategory][0]->Fill(trk_adc[iSpec][0][i]);
                h_hit_adc[iCategory][1]->Fill(trk_adc[iSpec][1][i]);
                h_hit_adc_ylocal[iCategory][0]->Fill(trk_y1_local[iSpec][i], trk_adc[iSpec][0][i]);
                h_hit_adc_ylocal[iCategory][1]->Fill(trk_y2_local[iSpec][i], trk_adc[iSpec][1][i]);
            }
            h_trk_has_hodo_hit->Fill(trk_has_hodo_hit[iSpec][i]);
            h_vertex_z[0]->Fill(trk_projz[iSpec][i], react_z[iSpec]);
            h_diff_x->Fill(trk_projx[iSpec][i] - react_x[iSpec]);
            h_diff_y->Fill(trk_projy[iSpec][i] - react_y[iSpec]);
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

}