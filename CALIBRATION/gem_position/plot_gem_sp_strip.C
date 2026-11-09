void plot_gem_sp_strip(){
     const int MAX_DATA = 10000; // Adjust as needed 
    TString inList = "filelist_23360_5sigma.txt";
    TString outName = "gem_sp_strip_23360_3sigma.root";

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
    TTreeReaderArray<Double_t>* istrip[nSpecs][2]; // [spec][gem]
    TTreeReaderArray<Double_t>* isU[nSpecs][2]; // [spec][gem]
    TTreeReaderArray<Double_t>* isV[nSpecs][2]; // [spec][gem]
    TTreeReaderArray<Double_t>* strip_adcsum[nSpecs][2]; // [spec][gem]
    TTreeReaderArray<Double_t>* strip_tmean[nSpecs][2]; // [spec][gem]
    TTreeReaderValue<Int_t>* nstrips[nSpecs][2]; // [spec][gem]
    TTreeReaderValue<Double_t>* react_x[nSpecs];
    TTreeReaderValue<Double_t>* react_y[nSpecs];
    TTreeReaderValue<Double_t>* react_z[nSpecs];
    TTreeReaderArray<Double_t>* sp_posX[nSpecs]; // [spec]
    TTreeReaderArray<Double_t>* sp_posY[nSpecs]; // [spec]
    TTreeReaderArray<Double_t>* sp_layer[nSpecs]; // [spec]
    TTreeReaderValue<Int_t>* nSp[nSpecs]; // [spec]

    for(int i = 0; i < nSpecs; i++) {
        for(int j = 0; j < 2; j++) { // Loop over the two GEMs
            istrip[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".gem.m%i.strip.istrip",j));
            isU[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".gem.m%i.strip.IsU",j));
            isV[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".gem.m%i.strip.IsV",j));
            strip_adcsum[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".gem.m%i.strip.ADCsum",j));
            nstrips[i][j] = new TTreeReaderValue<Int_t>(reader, Form("Ndata.%s.gem.m%i.strip.istrip", spec[i].Data(), j));
            strip_tmean[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".gem.m%i.strip.Tmean",j));
        }
        react_x[i] = new TTreeReaderValue<Double_t>(reader, spec[i] + ".react.x");
        react_y[i] = new TTreeReaderValue<Double_t>(reader, spec[i] + ".react.y");
        react_z[i] = new TTreeReaderValue<Double_t>(reader, spec[i] + ".react.z");
        sp_posX[i] = new TTreeReaderArray<Double_t>(reader, spec[i] + ".gem.sp.posX");
        sp_posY[i] = new TTreeReaderArray<Double_t>(reader, spec[i] + ".gem.sp.posY");
        sp_layer[i] = new TTreeReaderArray<Double_t>(reader, spec[i] + ".gem.sp.layer");
        nSp[i] = new TTreeReaderValue<Int_t>(reader, Form("Ndata.%s.gem.sp.posX", spec[i].Data()));
    }
     // Create histograms
    TH1D* h_strip[nSpecs][2][2]; // [spec][gem][U/V]
    TH2D* h_strip_adc_sum[nSpecs][2][2]; // [spec][gem][U/V]
    TH2D* h_strip_tmean[nSpecs][2][2]; // [spec][gem][U/V]
    TH1D* h_sp_posX[nSpecs][2]; // [spec][gem]
    TH1D* h_sp_posX_w[nSpecs][2]; // [spec][gem]
    TH1D* h_sp_posY[nSpecs][2]; // [spec][gem]
    TH1D* h_sp_posY_w[nSpecs][2]; // [spec][gem]
    TH1D* h_react_z[nSpecs]; // [spec]

    for(int i = 0; i < nSpecs; i++) {
        for(int j = 0; j < 2; j++) { // Loop over the two GEMs
            h_strip[i][j][0] = new TH1D(Form("h_strip_%s_gem%d_U", spec[i].Data(), j), Form("%s GEM%d U Strip Distribution", spec[i].Data(), j), 1536, -0.5, 1536.5);
            h_strip[i][j][1] = new TH1D(Form("h_strip_%s_gem%d_V", spec[i].Data(), j), Form("%s GEM%d V Strip Distribution", spec[i].Data(), j), 3072, -0.5, 3072.5);
            h_strip_adc_sum[i][j][0] = new TH2D(Form("h_strip_adc_sum_%s_gem%d_U", spec[i].Data(), j), Form("%s GEM%d U Strip ADC Sum Distribution", spec[i].Data(), j), 1536, -0.5, 1536.5, 1000, 0, 5000);
            h_strip_adc_sum[i][j][1] = new TH2D(Form("h_strip_adc_sum_%s_gem%d_V", spec[i].Data(), j), Form("%s GEM%d V Strip ADC Sum Distribution", spec[i].Data(), j), 3072, -0.5, 3072.5, 1000, 0, 5000);
            h_strip_tmean[i][j][0] = new TH2D(Form("h_strip_tmean_%s_gem%d_U", spec[i].Data(), j), Form("%s GEM%d U Strip Tmean Distribution", spec[i].Data(), j), 1536, -0.5, 1536.5, 130, 0, 130);
            h_strip_tmean[i][j][1] = new TH2D(Form("h_strip_tmean_%s_gem%d_V", spec[i].Data(), j), Form("%s GEM%d V Strip Tmean Distribution", spec[i].Data(), j), 3072, -0.5, 3072.5, 130, 0, 130);

            h_sp_posX[i][j] = new TH1D(Form("h_sp_posX_%s_gem%d", spec[i].Data(), j), Form("%s GEM%d SP PosX Distribution", spec[i].Data(), j), 100, -50, 50);
            h_sp_posY[i][j] = new TH1D(Form("h_sp_posY_%s_gem%d", spec[i].Data(), j), Form("%s GEM%d SP PosY Distribution", spec[i].Data(), j), 100, -50, 50);
            h_sp_posX_w[i][j] = new TH1D(Form("h_sp_posX_w_%s_gem%d", spec[i].Data(), j), Form("%s GEM%d SP PosX Weighted Distribution", spec[i].Data(), j), 100, -50, 50);
            h_sp_posY_w[i][j] = new TH1D(Form("h_sp_posY_w_%s_gem%d", spec[i].Data(), j), Form("%s GEM%d SP PosY Weighted Distribution", spec[i].Data(), j), 100, -50, 50);
        }
        h_react_z[i] = new TH1D(Form("h_react_z_%s", spec[i].Data()), Form("%s Reaction Z Distribution", spec[i].Data()), 100, -20, 20);
    }
     // Loop over entries
    int j = 0;
    while (reader.Next()) {
        int nHits[2]=  {0,0}; // [GEM0, GEM1]
       if (j++ % 1000 == 0) {
            std::cout << "Processed " << j << " entries.\r" << std::flush;
            /*if (j >= 1e6) {
                std::cout << "Reached 100k entries, stopping to save time." << std::endl;
                break; // Stop after processing 100k entries for testing
            }*/
        }
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
        h_react_z[iSpec]->Fill(**react_z[iSpec]);
        for(int j = 0; j < 2; j++) { // Loop over the two GEMs
            int nStrips = **nstrips[iSpec][j];
            for(int k = 0; k < nStrips; k++) {
                if ((*isU[iSpec][j])[k] == 1) {
                    h_strip[iSpec][j][0]->Fill((*istrip[iSpec][j])[k]);
                    h_strip_adc_sum[iSpec][j][0]->Fill((*istrip[iSpec][j])[k], (*strip_adcsum[iSpec][j])[k]);
                    h_strip_tmean[iSpec][j][0]->Fill((*istrip[iSpec][j])[k], (*strip_tmean[iSpec][j])[k]);
                }
                if ((*isV[iSpec][j])[k] == 1) {
                    h_strip[iSpec][j][1]->Fill((*istrip[iSpec][j])[k]);
                    h_strip_adc_sum[iSpec][j][1]->Fill((*istrip[iSpec][j])[k], (*strip_adcsum[iSpec][j])[k]);
                    h_strip_tmean[iSpec][j][1]->Fill((*istrip[iSpec][j])[k], (*strip_tmean[iSpec][j])[k]);
                }
            }
        }
        int nSp_event = **nSp[iSpec];
        for(int k = 0; k < nSp_event; k++) {
            if ((*sp_layer[iSpec])[k] == 0) { // Layer 0
                nHits[0]++;
                h_sp_posX[iSpec][0]->Fill((*sp_posX[iSpec])[k]);
                h_sp_posY[iSpec][0]->Fill((*sp_posY[iSpec])[k]);
            } else if ((*sp_layer[iSpec])[k] == 1) { // Layer 1
                nHits[1]++;
                h_sp_posX[iSpec][1]->Fill((*sp_posX[iSpec])[k]);
                h_sp_posY[iSpec][1]->Fill((*sp_posY[iSpec])[k]);
            }
        }
        for(int k=0; k<nSp_event; k++){
            if ((*sp_layer[iSpec])[k] == 0) { // Layer 0
                h_sp_posX_w[iSpec][0]->Fill((*sp_posX[iSpec])[k], nHits[1]);
                h_sp_posY_w[iSpec][0]->Fill((*sp_posY[iSpec])[k], nHits[1]);
            } else if ((*sp_layer[iSpec])[k] == 1) { // Layer 1
                h_sp_posX_w[iSpec][1]->Fill((*sp_posX[iSpec])[k], nHits[0]);
                h_sp_posY_w[iSpec][1]->Fill((*sp_posY[iSpec])[k], nHits[0]);
            }
        }
        
    }
        // Save histograms
    outputFile->cd();
    for(int i = 0; i < nSpecs; i++) {
        for(int j = 0; j < 2; j++) { // Loop over the two GEMs
            h_strip[i][j][0]->Write();
            h_strip[i][j][1]->Write();
            h_strip_adc_sum[i][j][0]->Write();
            h_strip_adc_sum[i][j][1]->Write();
            h_strip_tmean[i][j][0]->Write();
            h_strip_tmean[i][j][1]->Write();
            h_sp_posX[i][j]->Write();
            h_sp_posY[i][j]->Write();
            h_sp_posX_w[i][j]->Write();
            h_sp_posY_w[i][j]->Write();
        }
        h_react_z[i]->Write();
    }
    outputFile->Close();
     // Clean up
     for(int i = 0; i < nSpecs; i++) {
        for(int j = 0; j < 2; j++) { // Loop over the two GEMs
            delete istrip[i][j];
            delete isU[i][j];
            delete isV[i][j];
            delete nstrips[i][j];
            delete strip_adcsum[i][j];
            delete strip_tmean[i][j];
        }
        delete sp_posX[i];
        delete sp_posY[i];
        delete sp_layer[i];
        delete nSp[i];
        delete react_x[i];
        delete react_y[i];
        delete react_z[i];
    }
    delete outputFile;
    delete T;
    
}
