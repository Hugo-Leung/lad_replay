const int MAX_DATA  = 10000;
struct GEMHit{
    Double_t posX;
    Double_t posY;
};
void plotPos(){
    TString inList = "filelist_optics1.txt";
    const int nSpecs = 2;
    TString spec[] = {"H","P"}; // Spectrometer name (H or P)
    TChain* T = new TChain("T");
    std::ifstream fileList(inList);
    if (!fileList.is_open()) {
        std::cerr << "Error: Unable to open file list: " << inList << std::endl;
        return;
    }
    std::string line;
    while (std::getline(fileList, line)) {
        if (!line.empty()) {
            T->Add(line.c_str());
        }
    }
    fileList.close();

    TString outName = "gemHit_positions.root";
    TFile* outFile = new TFile(outName, "RECREATE"); 

    Long64_t nEntries = T->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;
    
    Double_t posX[2][MAX_DATA];
    Double_t posY[2][MAX_DATA];
    Double_t layer[2][MAX_DATA];
    Int_t nPointsX[2];
    Int_t nPointsY[2];
    Int_t nLayers[2];
    UInt_t triggerBits=0;
    Double_t triggerBits_float=0;
    UInt_t triggerMask[]={8,2};// HMS Ereal, SHMS Ereal

    // Set branch addresses for HMS and SHMS
    for (int iSpec=0; iSpec<nSpecs; iSpec++){
        T->SetBranchAddress(spec[iSpec]+".gem.sp.posX", &posX[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.sp.posY", &posY[iSpec]);
        T->SetBranchAddress("Ndata."+spec[iSpec]+".gem.sp.posX", &nPointsX[iSpec]);
        T->SetBranchAddress("Ndata."+spec[iSpec]+".gem.sp.posY", &nPointsY[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.sp.layer", &layer[iSpec]);
        T->SetBranchAddress("Ndata."+spec[iSpec]+".gem.sp.layer", &nLayers[iSpec]);
    }

    TH2D* h2_gemPos_layer0 = new TH2D("h2_gemPos_layer0", "GEM Hit Positions Layer 0;X Position (cm);Y Position (cm)", 100, 20, 120, 100, -35, 35);
    TH2D* h2_gemPos_layer1 = new TH2D("h2_gemPos_layer1", "GEM Hit Positions Layer 1;X Position (cm);Y Position (cm)", 100, 20, 120, 100, -35, 35);
    TH1D* h_locationCut = new TH1D("h_locationCut", "Location Cut Status;Status;Counts", 2, -0.5, 1.5);
    TH2D* h2_locationCut_layer0 = new TH2D("h2_locationCut_layer0", "GEM Hit Positions Layer 0 with Location Cut;X Position (cm);Y Position (cm)", 100, 20, 120, 100, -35, 35);
    TH2D* h2_locationCut_layer1 = new TH2D("h2_locationCut_layer1", "GEM Hit Positions Layer 1 with Location Cut;X Position (cm);Y Position (cm)", 100, 20, 120, 100, -35, 35);
    T->SetBranchAddress("g.trigbits", &triggerBits_float);
    std::vector<GEMHit> gemHitsLayer0;
    std::vector<GEMHit> gemHitsLayer1;
    for (int i=0; i<nEntries; i++){
        T->GetEntry(i);
        int iSpec=-1;
        triggerBits = static_cast<UInt_t>(triggerBits_float);
        if ((triggerBits&triggerMask[0])!=0){//HMS
            iSpec=0;
        }
        if ((triggerBits&triggerMask[1])!=0){//SHMS
            iSpec=1;
        }
        if (iSpec==-1){
            continue; // no valid spectrometer triggered
        }
        if(nPointsX[iSpec]!=nPointsY[iSpec]||nPointsX[iSpec]!=nLayers[iSpec]) {
            cout<<"inconsistent number of points"<<std::endl;
            continue;
        }
        
        //if (i>=1e5) break; // limit number of entries to process for testing
        bool locationCut = false;
        for (int j=0; j<nPointsX[iSpec]; j++){
            //std::cout<<"Spec: "<<spec[iSpec]<<" Entry: "<<i<<" Point: "<<j<<" Layer: "<<layer[iSpec][j]<<" PosX: "<<posX[iSpec][j]<<" PosY: "<<posY[iSpec][j]<<std::endl;
            if (layer[iSpec][j]==0){
                GEMHit hit;
                hit.posX = posX[iSpec][j];
                hit.posY = posY[iSpec][j];
                if (posX[iSpec][j]>70 && posX[iSpec][j]<75 && posY[iSpec][j]>-2.5 && posY[iSpec][j]<2.5){
                    locationCut = true;
                    gemHitsLayer0.push_back(hit);
                }
                h2_gemPos_layer0->Fill(hit.posX, hit.posY);
            }
            else if (layer[iSpec][j]==1){
                GEMHit hit;
                hit.posX = posX[iSpec][j];
                hit.posY = posY[iSpec][j];
                gemHitsLayer1.push_back(hit);
            }
        }
        if(gemHitsLayer0.size()<2){
            locationCut = false;
        }
        for (const auto& hit0 : gemHitsLayer0) {
            h2_locationCut_layer0->Fill(hit0.posX, hit0.posY);
        }
        for (const auto& hit1 : gemHitsLayer1) {
            if (locationCut&&gemHitsLayer1.size()<2){
                h2_locationCut_layer1->Fill(hit1.posX, hit1.posY);
            }
            h2_gemPos_layer1->Fill(hit1.posX, hit1.posY); 
        }
        if (locationCut){
            h_locationCut->Fill(1);
        }else{
            h_locationCut->Fill(0);
        }
        gemHitsLayer0.clear();
        gemHitsLayer1.clear();
        if (i%1000==0){
            std::cout<<"\r processed "<<i<<"/"<<nEntries<<std::flush;
        }
    }
    std::cout<< " finished processing all entries"<<std::endl;
    outFile->cd();
    h2_gemPos_layer0->Write();
    h2_gemPos_layer1->Write();
    h_locationCut->Write();
    h2_locationCut_layer0->Write();
    h2_locationCut_layer1->Write();
    TH2D* h_ratio = (TH2D*)h2_locationCut_layer1->Clone("h_ratio_locationCut_layer1");;
    h_ratio->Divide(h2_gemPos_layer1);
    h_ratio->SetTitle("Ratio of GEM Layer 1 Hits after Location Cut to All Hits;X Position (cm);Y Position (cm)");
    h_ratio->Write();
    outFile->Close();
    delete T;



    

}