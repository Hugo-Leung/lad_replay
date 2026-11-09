const double plane_r[5]     = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935}; // Radius of the planes in cm
const double c  = 29.9792458; // speed of light in cm/ns

void plot_lad_tof_simple(){
    TString inList = "filelist_23367_23370.txt";
    TString outName = "lad_tof_simple_suppressfirstTS.root";
    TChain *C = new TChain("T");
    ifstream inFile(inList);
    TString currentLine;
    while(currentLine.ReadLine(inFile)){
        if(!currentLine.BeginsWith("#")){
            C->Add(currentLine);
        }
    }

    TFile *outputFile = new TFile(outName, "RECREATE");

	TString spec[] = {"H","P"}; // Spectrometer name (H or P)
	UInt_t triggerMask[]={8,2};// HMS Ereal, SHMS Ereal
	const int nSpecs = 2;
    TTreeReader reader(C);
    TTreeReaderValue<Double_t>* triggerBits_float = new TTreeReaderValue<Double_t>(reader, "g.trigbits");
    TTreeReaderArray<Double_t>* tof_rfcorr[nSpecs][2]; // front and back TOF
    TTreeReaderValue<Double_t>* nGoodhit[nSpecs]; // number of TOF hits
    TTreeReaderArray<Double_t>* isProton[nSpecs][2]; // is proton flag
    TTreeReaderArray<Double_t>* paddle[nSpecs][2]; // paddle number
    TTreeReaderArray<Double_t>* plane[nSpecs][2]; // plane number
    TTreeReaderArray<Double_t>* trackID[nSpecs];// track ID
    TTreeReaderArray<Double_t>* ypos[nSpecs][2]; // y position of the hit

    for (int i = 0; i < nSpecs; i++) {
        for (int j = 0; j < 2; j++) { // Loop over the two hodo planes
            tof_rfcorr[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".ladhod.goodhit_hit_tof_%d", j));
            isProton[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".ladhod.goodhit_isProton_%d", j));
            paddle[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".ladhod.goodhit_paddle_%d", j));
            plane[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".ladhod.goodhit_plane_%d", j));
            ypos[i][j] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".ladhod.goodhit_hit_ypos_%d", j));
        }
        nGoodhit[i] = new TTreeReaderValue<Double_t>(reader, spec[i] + Form(".ladhod.goodhit_n"));
        trackID[i] = new TTreeReaderArray<Double_t>(reader, spec[i] + Form(".ladhod.goodhit_trackid"));
    }
    double tofRange[2]={-150,175}; // ns
    int nBins=650;
    TH1F *h_tof_rfcorr[nSpecs][2][3]; // front and back TOF, all hits, with trackID, without trackID
    for (int i = 0; i < nSpecs; i++) {
        for (int j = 0; j < 2; j++) {
            h_tof_rfcorr[i][j][0] = new TH1F(Form("h_tof_rfcorr_%s_plane%d_all", spec[i].Data(), j), Form("TOF RF-corrected %s plane %d (all hits);TOF [ns];Counts", spec[i].Data(), j), nBins, tofRange[0], tofRange[1]);
            h_tof_rfcorr[i][j][1] = new TH1F(Form("h_tof_rfcorr_%s_plane%d_withTrackID", spec[i].Data(), j), Form("TOF RF-corrected %s plane %d (with trackID);TOF [ns];Counts", spec[i].Data(), j), nBins, tofRange[0], tofRange[1]);
            h_tof_rfcorr[i][j][2] = new TH1F(Form("h_tof_rfcorr_%s_plane%d_withoutTrackID", spec[i].Data(), j), Form("TOF RF-corrected %s plane %d (without trackID);TOF [ns];Counts", spec[i].Data(), j), nBins, tofRange[0], tofRange[1]);
        }
    }
    int count = 0;
    while (reader.Next()) {
        if (count++ % 10000 == 0) {
            std::cout << "\rProcessing " << count << " events." << std::flush;
        }
        /*if (count > 1e6) {
            std::cout << "\nReached 1 million events, stopping to save time." << std::endl;
            break; // Stop after processing 1 million events for testing
        }*/
        UInt_t triggerBits = static_cast<UInt_t>(**triggerBits_float);
        int iSpec = -1;
        if ((triggerBits & triggerMask[0]) != 0) { // HMS
            iSpec = 0;
        } else if ((triggerBits & triggerMask[1]) != 0) { // SHMS
            iSpec = 1;
        } else {
            std::cerr << "Warning: Event does not match expected trigger conditions. Skipping event." << std::endl;
            std::cout << "Trigger bits: " << (triggerBits & 0xff) << std::endl;
            continue; // Skip events that don't match the trigger conditions
        }
        int nHits = static_cast<int>(**nGoodhit[iSpec]);
        for (int k = 0; k < nHits; k++) {
            double tof = (*tof_rfcorr[iSpec][0])[k];
            if (tof>tofRange[0]-2000 && tof<tofRange[1]+2000 && (*plane[iSpec][0])[k] >= 0 ) { // front plane
                double pathlength = plane_r[(int)(*plane[iSpec][0])[k]]; // cm
                pathlength= pathlength*pathlength + (*ypos[iSpec][0])[k]*(*ypos[iSpec][0])[k];
                pathlength += (110-22*(*paddle[iSpec][0])[k])*(110-22*(*paddle[iSpec][0])[k]); // distance from target to front plane
                pathlength = sqrt(pathlength);
                double lightTime = pathlength / c; // ns
                double correctedTOF = tof - lightTime;
                //h_tof_rfcorr[iSpec][0][0]->Fill(tof); // all hits
                h_tof_rfcorr[iSpec][0][0]->Fill(correctedTOF); // all hits
                if ((*trackID[iSpec])[k] >= 0) {
                    h_tof_rfcorr[iSpec][0][1]->Fill(correctedTOF); // with trackID
                } else {
                    h_tof_rfcorr[iSpec][0][2]->Fill(correctedTOF); // without trackID
                }
            }
            tof = (*tof_rfcorr[iSpec][1])[k];
            if (tof>tofRange[0]-2000 && tof<tofRange[1]+2000 && (*plane[iSpec][1])[k] >= 0 ) { // back plane
                double pathlength = plane_r[(int)(*plane[iSpec][1])[k]]; // cm
                pathlength= pathlength*pathlength + (*ypos[iSpec][1])[k]*(*ypos[iSpec][1])[k];
                pathlength += (110-22*(*paddle[iSpec][1])[k])*(110-22*(*paddle[iSpec][1])[k]); // distance from target to back plane
                pathlength = sqrt(pathlength);
                double lightTime = pathlength / c; // ns
                double correctedTOF = tof - lightTime;
                //h_tof_rfcorr[iSpec][1][0]->Fill(tof); // all hits
                h_tof_rfcorr[iSpec][1][0]->Fill(correctedTOF); // all hits
                if ((*trackID[iSpec])[k] >= 0) {
                    h_tof_rfcorr[iSpec][1][1]->Fill(correctedTOF); // with trackID
                } else {
                    h_tof_rfcorr[iSpec][1][2]->Fill(correctedTOF); // without trackID
                }
            }
        }
    }

    // Save histograms
    outputFile->cd();
    for (int i = 0; i < nSpecs; i++) {
        for (int j = 0; j < 2; j++) { // Loop over the two hodo planes
            for (int k = 0; k < 3; k++) { // Loop over all hits, with trackID, without trackID
                h_tof_rfcorr[i][j][k]->Write();
            }
        }
    }
    //clean up and close the output file
    outputFile->Close();
    delete triggerBits_float;
    for (int i = 0; i < nSpecs; i++) {
        for (int j = 0; j < 2; j++) {
            delete tof_rfcorr[i][j];
            delete isProton[i][j];
            delete paddle[i][j];
            delete plane[i][j];
        }
        delete nGoodhit[i];
        delete trackID[i];
    }
    delete C;
}