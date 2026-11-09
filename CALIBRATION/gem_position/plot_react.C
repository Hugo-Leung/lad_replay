void plot_react(){
    TString inList = "filelist_3pass_optics1.txt";
    TString outName = "react_3pass_optics1.root";

    /*double fixedz[]={-5,5};//optics2;
    double z_range[] = {-10,0,10};
    const int nFixedz =2; */
    double fixedz[] = {-10, 0, 10};
    double z_range[] = {-15, -5, 5, 15}; // Define the z ranges for the fixed positions, should be one more than the number of fixed z positions
    const int nFixedz = 3;

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

    Double_t react_x[2], react_y[2], react_z[2];
    Double_t dp[2], etracknorm[2], npeSum[2];

    UInt_t triggerBits=0;
    Double_t triggerBits_float=0;
    UInt_t triggerMask[]={8,2};// HMS Ereal, SHMS Ereal
    T->SetBranchAddress("g.trigbits", &triggerBits_float);
    for (int iSpec=0; iSpec<nSpecs; iSpec++){
        T->SetBranchAddress(spec[iSpec]+".react.x", &react_x[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".react.y", &react_y[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".react.z", &react_z[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gtr.dp", &dp[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".cal.etracknorm", &etracknorm[iSpec]);
        if(iSpec==0){
           T->SetBranchAddress(spec[iSpec]+".cer.npeSum", &npeSum[iSpec]); 
        }else if (iSpec==1){
           T->SetBranchAddress(spec[iSpec]+".ngcer.npeSum", &npeSum[iSpec]);  
        }
    }
    const int nHisto = 2;
    TH1D* h_vertex_z[nHisto];

    int nEntries = T->GetEntries();
    for(int i=0; i<nHisto; i++){
       h_vertex_z[i]=new TH1D(Form("h_vertex_z_%d", i), ";"+spec[i]+".react.z",100,-15,15); 
    }
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
            continue; // no valid spectrometer triggered
        }

        if (iSpec==0){
            if(TMath::Abs(dp[iSpec])>=10 || etracknorm[iSpec]<=0.65 || npeSum[iSpec]<=6)
            continue;
        }else if (iSpec==1){
            if(dp[iSpec]>20 || dp[iSpec]<-15 || etracknorm[iSpec]<=0.8 || npeSum[iSpec]<=5)
            continue;
        }
        else {
            continue;
        }
        h_vertex_z[iSpec]->Fill(react_z[iSpec]);   
    }
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
            TFitResultPtr fitResult = h_vertex_z[i]->Fit("gaus", "QLISNR", "", fixedz[ifoil]-2.5, fixedz[ifoil]+2.5);
            double mean = fitResult->Parameter(1);
            double sigma = fitResult->Parameter(2);
            cout<<"Z = "<<fixedz[ifoil]<<": Mean = "<<mean<<" cm, Sigma = "<<sigma<<" cm"<<std::endl;
            double fitRangeMin = TMath::Max (mean - 1.5*sigma,-25.0);
            double fitRangeMax = TMath::Min (mean + 1.5*sigma, 25.0); 
            TFitResultPtr fitResult2 = h_vertex_z[i]->Fit("gaus", "QLIS+R", "", fitRangeMin, fitRangeMax);
            means[ifoil] = fitResult2->Parameter(1);
            means_err[ifoil] = fitResult2->ParError(1);
            sigmas[ifoil] = fitResult2->Parameter(2);
            sigmas_err[ifoil] = fitResult2->ParError(2);
            pt->AddText(Form("Z = %.1f cm: Mean = %.2f #pm %.2f cm, Sigma = %.2f #pm %.2f cm", fixedz[ifoil], means[ifoil], means_err[ifoil], sigmas[ifoil], sigmas_err[ifoil]));
        }
        h_vertex_z[i]->GetListOfFunctions()->Add(pt);
    
    }
    outputFile->Write();
    outputFile->Close();
}