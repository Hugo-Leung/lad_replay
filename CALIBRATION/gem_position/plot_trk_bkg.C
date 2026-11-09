void plot_trk_bkg(){
    //get a Out-of-time removed vertex distribution from LAD

    
    TString inList = "filelist_23367_23370.txt";
    TString outName = "gem_trk_bkgSub_23367_23370_tCut_noZS.root";
    TString spec[] = {"H","P"}; // Spectrometer name (H or P)
    const int nSpecs = 2;

    const double hodoTime[2] = {1600, 2010}; //Time window for hodo hit
    const double hodoTimeSignal[2] = {1670, 1875}; //Time window for signal hodo hit
    const double chisqCut[2] = {20, 25}; //Chi2 cut for different categories
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
    TTreeReaderValue<Double_t>* nGEM2Dhits[2][2];
    TTreeReaderArray<Double_t>* trk_spid[2][2];
    TTreeReaderArray<Double_t>* sp_layer[2];
    TTreeReaderArray<Double_t>* sp_time[2];
    TTreeReaderArray<Double_t>* sp_spID[2];
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
    TTreeReaderArray<Double_t>* hitTof_rfcorr[2][2];
    TTreeReaderArray<Double_t>* hitTrkID[2][2];
    TTreeReaderArray<Double_t>* trk_adc[2][2];
    for (int i=0; i<2; i++){
        nGEM2Dhits[i][0] = new TTreeReaderValue<Double_t>(reader, Form("%s.gem.m%i.n2Dhits", spec[i].Data(), 0));
        nGEM2Dhits[i][1] = new TTreeReaderValue<Double_t>(reader, Form("%s.gem.m%i.n2Dhits", spec[i].Data(), 1));
        trk_spid[i][0] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.spID_0", spec[i].Data()));
        trk_spid[i][1] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.trk.spID_1", spec[i].Data()));

        sp_layer[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.sp.layer", spec[i].Data()));
        sp_spID[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.sp.spID", spec[i].Data()));
        sp_time[i] = new TTreeReaderArray<Double_t>(reader, Form("%s.gem.sp.time", spec[i].Data()));
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
        hitTof_rfcorr[i][0] = new TTreeReaderArray<Double_t>(reader, Form("%s.ladhod.goodhit_hit_tof_rfcorr_0", spec[i].Data()));
        hitTof_rfcorr[i][1] = new TTreeReaderArray<Double_t>(reader, Form("%s.ladhod.goodhit_hit_tof_rfcorr_1", spec[i].Data()));
    }
    

    //Define histograms
    //categories: 0: all good tracks, 1: no hodo hit, 2: 1 hodo hit, 3: 2 hodo hits
    //spectrometers: 0: H, 1: P
    // intime, out-of-time, and background distributions for each spectrometer
    const int nHisto = 2; //1 hodo hit, 2 hodo hits, don't really care about no hodo hit for this analysis
    const int nTimePeriods = 4; // all, in-time, out-of-time, subtracted
    TString timePeriods[4] = {"all", "inTime", "outOfTime", "subtracted"};
    TString hodoNames[2] = { "oneHodoHit", "twoHodoHits"};
    TH1D* h_trk_projz[nHisto][nSpecs][nTimePeriods];
    TH1D* h_trk_chisq[nHisto][nSpecs][nTimePeriods];
    TH1D* h_trk_theta[nHisto][nSpecs][nTimePeriods];
    TH1D* h_trk_phi[nHisto][nSpecs][nTimePeriods];
    TH1D* h_hodoTime[nHisto][nSpecs][nTimePeriods][2]; // 2 planes of hodoscope
    TH1D* h_hodoTof[nHisto][nSpecs][nTimePeriods][2]; // 2 planes of hodoscope
    TH1D* h_gemTime[nHisto][nSpecs][nTimePeriods][2]; // 2 planes of GEM

    for (int iCategory=0; iCategory<nHisto; iCategory++){
        for (int iSpec=0; iSpec<nSpecs; iSpec++){
            for (int iTime=0; iTime<nTimePeriods; iTime++){
                h_trk_projz[iCategory][iSpec][iTime] = new TH1D(Form("h_trk_projz_%s_%s_%s", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), Form("Track projected z distribution (%s, %s, %s)", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), 100, -20, 20);
                h_trk_chisq[iCategory][iSpec][iTime] = new TH1D(Form("h_trk_chisq_%s_%s_%s", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), Form("Track chi-squared distribution (%s, %s, %s)", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), 100, 0, 20);
                h_trk_theta[iCategory][iSpec][iTime] = new TH1D(Form("h_trk_theta_%s_%s_%s", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), Form("Track theta distribution (%s, %s, %s)", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), 100, 70, 170);
                h_trk_phi[iCategory][iSpec][iTime] = new TH1D(Form("h_trk_phi_%s_%s_%s", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), Form("Track phi distribution (%s, %s, %s)", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), 100, -35, 35);
                
                h_gemTime[iCategory][iSpec][iTime][0] = new TH1D(Form("h_gemTime_plane0_%s_%s_%s", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), Form("GEM plane 0 hit time distribution (%s, %s, %s)", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), 140, 0, 140);
                h_gemTime[iCategory][iSpec][iTime][1] = new TH1D(Form("h_gemTime_plane1_%s_%s_%s", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), Form("GEM plane 1 hit time distribution (%s, %s, %s)", hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), 140, 0, 140);


                //sumw2
                h_trk_projz[iCategory][iSpec][iTime]->Sumw2();
                h_trk_chisq[iCategory][iSpec][iTime]->Sumw2();
                h_trk_theta[iCategory][iSpec][iTime]->Sumw2();
                h_trk_phi[iCategory][iSpec][iTime]->Sumw2();

                h_gemTime[iCategory][iSpec][iTime][0]->Sumw2();
                h_gemTime[iCategory][iSpec][iTime][1]->Sumw2();

                for (int plane=0; plane<2; plane++){
                    h_hodoTime[iCategory][iSpec][iTime][plane] = new TH1D(Form("h_hodoTime_plane%d_%s_%s_%s", plane, hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), Form("Hodoscope hit time distribution (plane %d, %s, %s, %s)", plane, hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), (hodoTime[1]-hodoTime[0])/10, hodoTime[0], hodoTime[1]);
                    h_hodoTime[iCategory][iSpec][iTime][plane]->Sumw2();
                    h_hodoTof[iCategory][iSpec][iTime][plane] = new TH1D(Form("h_hodoTof_plane%d_%s_%s_%s", plane, hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), Form("Hodoscope hit time-of-flight distribution (plane %d, %s, %s, %s)", plane, hodoNames[iCategory].Data(), spec[iSpec].Data(), timePeriods[iTime].Data()), 1000, -100, 300);
                    h_hodoTof[iCategory][iSpec][iTime][plane]->Sumw2();
                }
            }
        }
    }
    // Loop over entries in the TTree and fill histograms
    int j = 0;
    while (reader.Next()) {

        if (j++ % 1000 == 0) {
            std::cout << "Processed " << j << " entries.\r" << std::flush;
        }
        //if (j>100000) break; // Limit to 500,000 entries for testing
        UInt_t triggerBits = static_cast<UInt_t>(**triggerBits_float);
        bool hmsEreal = (triggerBits & 8) != 0;
        bool shmsEreal = (triggerBits & 2) != 0;
        for (int iSpec=0; iSpec<nSpecs; iSpec++){
           
            int nTracks_iSpec = **nTrack[iSpec];
            for (int i = 0; i < nTracks_iSpec; i++) {
                int iCategory = -1; // Default category for all good tracks
                if ((*trk_has_hodo_hit[iSpec])[i] == 1) {
                    iCategory = 0; // 1 hodo hit
                } else if ((*trk_has_hodo_hit[iSpec])[i] == 2) {
                    iCategory = 1; // 2  hodo hits
                }else {
                    continue; // Skip tracks that don't fall into the defined categories
                }
                Double_t hitTime_tmp[2] = {-1, -1}; // Temporary array to store hit times for the current track
                Double_t hitTof_tmp[2] = {-1, -1}; // Temporary array to store hit times for the current track
                Double_t gemTime_tmp[2] = {-1, -1}; // Temporary array to store GEM times for the current track
                for (int j=0; j<**nGoodHits[iSpec]; j++){
                    for (int plane=0; plane<2; plane++){
                        if ((*hitTrkID[iSpec][0])[j]==i){
                            hitTime_tmp[plane]=(*hitTime[iSpec][plane])[j];
                            hitTof_tmp[plane]=(*hitTof_rfcorr[iSpec][plane])[j];
                        }
                    }
                }
                //layer 0
                int spIDX= (*trk_spid[iSpec][0])[i];
                if (spIDX>=0 && spIDX<**nGEM2Dhits[iSpec][0]&& (*sp_layer[iSpec])[spIDX]==0 && (*sp_spID[iSpec])[spIDX]-(*trk_spid[iSpec][0])[i]<0.1){
                    gemTime_tmp[0] = (*sp_time[iSpec])[spIDX];
                }else{
                    cout<< "Warning: Track " << i << " in spectrometer " << spec[iSpec].Data() << " has invalid GEM layer 0 index: " << spIDX << std::endl;
                }
                //layer 1
                spIDX= (*trk_spid[iSpec][1])[i]+ (**nGEM2Dhits[iSpec][0]);
                if (spIDX>=(**nGEM2Dhits[iSpec][0]) && spIDX<**nGEM2Dhits[iSpec][1]+(**nGEM2Dhits[iSpec][0]) && fabs((*sp_layer[iSpec])[spIDX]-1)<0.01 && (*sp_spID[iSpec])[spIDX]-(*trk_spid[iSpec][1])[i]<0.1){
                    gemTime_tmp[1] = (*sp_time[iSpec])[spIDX];
                }else{
                    cout<<"Track "<< i << "in Spec" << spec[iSpec].Data() <<": SpID "<< (*sp_spID[iSpec])[spIDX] << " != track_spID1 "<< (*trk_spid[iSpec][1])[i]<< " or layerID: "<< (*sp_layer[iSpec])[spIDX] << std::endl;
                }


                if ((*trk_is_good[iSpec])[i] != 1) {
                    continue; // Skip tracks that are not marked as good
                }
                //chisq cut and z cut for all categories
                if(!((iCategory==0&&(*trk_chisq[iSpec])[i]<chisqCut[0]) || (iCategory==1&&(*trk_chisq[iSpec])[i]<chisqCut[1]))||((*trk_projz[iSpec])[i]<-19 || (*trk_projz[iSpec])[i]>19)){
                    continue; // Skip tracks that don't meet the chisq and z cut criteria
                }
                int iTimePeriod = 0; // Default to all tracks
                if (iCategory==0){
                    //only 1 hodo hit, check if in-time or out-of-time
                    bool accept=false;
                    if (hitTime_tmp[0]>hodoTimeSignal[0] && hitTime_tmp[0]<hodoTimeSignal[1]){
                        iTimePeriod = 1; //in-time
                        accept=true;
                    } else if (hitTime_tmp[0]>hodoTime[0] && hitTime_tmp[0]<hodoTime[1]){
                        iTimePeriod = 2; //out-of-time
                        accept=true;
                    } 
                    if (!accept&&hitTime_tmp[1]>0){
                        //check if the second plane has a hit
                        if (hitTime_tmp[1]>hodoTimeSignal[0] && hitTime_tmp[1]<hodoTimeSignal[1]){
                            iTimePeriod = 1; //in-time
                            accept=true;
                        } else if (hitTime_tmp[1]>hodoTime[0] && hitTime_tmp[1]<hodoTime[1]){
                            iTimePeriod = 2; //out-of-time
                            accept=true;
                        }

                    }
                    if (!accept){
                        continue; //skip tracks that are not in-time or out-of-time
                    }
                } else if (iCategory==1){
                    bool accept=false;
                    //require both planes to be in-time or out-of-time
                    if ((hitTime_tmp[0]>hodoTimeSignal[0] && hitTime_tmp[0]<hodoTimeSignal[1]) && (hitTime_tmp[1]>hodoTimeSignal[0] && hitTime_tmp[1]<hodoTimeSignal[1])){
                        iTimePeriod = 1; //in-time
                        accept=true;
                    } else if ((hitTime_tmp[0]>hodoTime[0] && hitTime_tmp[0]<hodoTime[1]) && (hitTime_tmp[1]>hodoTime[0] && hitTime_tmp[1]<hodoTime[1])){
                        iTimePeriod = 2; //out-of-time
                        accept=true;
                    }
                    if (!accept){
                        continue; //skip tracks that are not in-time or out-of-time
                    }
                }


                //fill histograms, both all (0) and the specific time period (iTimePeriod)
                h_trk_projz[iCategory][iSpec][0]->Fill((*trk_projz[iSpec])[i]);
                h_trk_projz[iCategory][iSpec][iTimePeriod]->Fill((*trk_projz[iSpec])[i]);
                h_trk_chisq[iCategory][iSpec][0]->Fill((*trk_chisq[iSpec])[i]);
                h_trk_chisq[iCategory][iSpec][iTimePeriod]->Fill((*trk_chisq[iSpec])[i]);
                h_trk_theta[iCategory][iSpec][0]->Fill((*trk_theta[iSpec])[i]*TMath::RadToDeg());
                h_trk_theta[iCategory][iSpec][iTimePeriod]->Fill((*trk_theta[iSpec])[i]*TMath::RadToDeg());
                h_trk_phi[iCategory][iSpec][0]->Fill((*trk_phi[iSpec])[i]*TMath::RadToDeg());
                h_trk_phi[iCategory][iSpec][iTimePeriod]->Fill((*trk_phi[iSpec])[i]*TMath::RadToDeg());
                for (int plane=0; plane<2; plane++){
                    if (hitTime_tmp[plane]>0){
                        h_hodoTime[iCategory][iSpec][0][plane]->Fill(hitTime_tmp[plane]);
                        h_hodoTime[iCategory][iSpec][iTimePeriod][plane]->Fill(hitTime_tmp[plane]);
                        h_hodoTof[iCategory][iSpec][0][plane]->Fill(hitTof_tmp[plane]);
                        h_hodoTof[iCategory][iSpec][iTimePeriod][plane]->Fill(hitTof_tmp[plane]);
                    }
                    h_gemTime[iCategory][iSpec][0][plane]->Fill(gemTime_tmp[plane]);
                    h_gemTime[iCategory][iSpec][iTimePeriod][plane]->Fill(gemTime_tmp[plane]);
                }

            }
        }

    }
    //calculate background subtracted distributions
    int iTimePeriod=3; //background subtracted
    const double scaleFactor = (hodoTimeSignal[1]-hodoTimeSignal[0])/((hodoTime[1]-hodoTime[0])-(hodoTimeSignal[1]-hodoTimeSignal[0]));
    cout << "Scale factor for background subtraction: " << scaleFactor << std::endl;
    for (int iCategory=0; iCategory<nHisto; iCategory++){
        for (int iSpec=0; iSpec<nSpecs; iSpec++){
            h_trk_projz[iCategory][iSpec][iTimePeriod]->Add(h_trk_projz[iCategory][iSpec][1], 1);
            h_trk_projz[iCategory][iSpec][iTimePeriod]->Add(h_trk_projz[iCategory][iSpec][2], -scaleFactor);

            h_trk_chisq[iCategory][iSpec][iTimePeriod]->Add(h_trk_chisq[iCategory][iSpec][1], 1);
            h_trk_chisq[iCategory][iSpec][iTimePeriod]->Add(h_trk_chisq[iCategory][iSpec][2], -scaleFactor);

            h_trk_theta[iCategory][iSpec][iTimePeriod]->Add(h_trk_theta[iCategory][iSpec][1], 1);
            h_trk_theta[iCategory][iSpec][iTimePeriod]->Add(h_trk_theta[iCategory][iSpec][2], -scaleFactor);

            h_trk_phi[iCategory][iSpec][iTimePeriod]->Add(h_trk_phi[iCategory][iSpec][1], 1);
            h_trk_phi[iCategory][iSpec][iTimePeriod]->Add(h_trk_phi[iCategory][iSpec][2], -scaleFactor);

            for (int plane=0; plane<2; plane++){
                h_hodoTime[iCategory][iSpec][iTimePeriod][plane]->Add(h_hodoTime[iCategory][iSpec][1][plane], 1);
                h_hodoTime[iCategory][iSpec][iTimePeriod][plane]->Add(h_hodoTime[iCategory][iSpec][2][plane], -scaleFactor);
                h_hodoTof[iCategory][iSpec][iTimePeriod][plane]->Add(h_hodoTof[iCategory][iSpec][1][plane], 1);
                h_hodoTof[iCategory][iSpec][iTimePeriod][plane]->Add(h_hodoTof[iCategory][iSpec][2][plane], -scaleFactor);
                h_gemTime[iCategory][iSpec][iTimePeriod][plane]->Add(h_gemTime[iCategory][iSpec][1][plane], 1);
                h_gemTime[iCategory][iSpec][iTimePeriod][plane]->Add(h_gemTime[iCategory][iSpec][2][plane], -scaleFactor);
            }
        }
    }

    //add vertical lines to hodo time histograms to indicate the in-time window
    for (int iCategory=0; iCategory<nHisto; iCategory++){
        for (int iSpec=0; iSpec<nSpecs; iSpec++){
            iTimePeriod=0; // only add lines to the "all" time period histograms
            for (int plane=0; plane<2; plane++){
                TLine* line1 = new TLine(hodoTimeSignal[0], 0, hodoTimeSignal[0], h_hodoTime[iCategory][iSpec][iTimePeriod][plane]->GetMaximum());
                line1->SetLineColor(kRed);
                line1->SetLineWidth(2);
                h_hodoTime[iCategory][iSpec][iTimePeriod][plane]->GetListOfFunctions()->Add(line1);
                TLine* line2 = new TLine(hodoTimeSignal[1], 0, hodoTimeSignal[1], h_hodoTime[iCategory][iSpec][iTimePeriod][plane]->GetMaximum());
                line2->SetLineColor(kRed);
                line2->SetLineWidth(2);
                h_hodoTime[iCategory][iSpec][iTimePeriod][plane]->GetListOfFunctions()->Add(line2);
            }
        }
    }

    //add horizontal line at y=0 to the background subtracted z distribution histograms
    for (int iCategory=0; iCategory<nHisto; iCategory++){
        for (int iSpec=0; iSpec<nSpecs; iSpec++){
            iTimePeriod=3; // only add line to the background subtracted histograms
            TLine* line = new TLine(h_trk_projz[iCategory][iSpec][iTimePeriod]->GetXaxis()->GetXmin(), 0, h_trk_projz[iCategory][iSpec][iTimePeriod]->GetXaxis()->GetXmax(), 0);
            line->SetLineColor(kBlack);
            line->SetLineWidth(1);
            h_trk_projz[iCategory][iSpec][iTimePeriod]->GetListOfFunctions()->Add(line);
            double error = 0;
            double integral = h_trk_projz[iCategory][iSpec][iTimePeriod]->IntegralAndError(1, h_trk_projz[iCategory][iSpec][iTimePeriod]->GetNbinsX()+1, error);
            cout << "Background subtracted integral for " << hodoNames[iCategory].Data() << " in " << spec[iSpec].Data() << ": " << integral << " +/- " << error << std::endl;
        }
    }
    //Write histograms to output file and clean up
    outputFile->Write();
    outputFile->Close();
    delete triggerBits_float;
    for (int i=0; i<2; i++){
        delete nTrack[i];
        delete nGEM2Dhits[i][0];
        delete nGEM2Dhits[i][1];
        delete trk_spid[i][0];
        delete trk_spid[i][1];
        delete sp_layer[i];
        delete sp_time[i];
        delete sp_spID[i];
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