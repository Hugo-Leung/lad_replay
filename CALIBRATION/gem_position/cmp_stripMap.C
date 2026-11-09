void cmp_stripMap(){
    TString inName[2] = {"filelist_optics1_v1.txt", "filelist_optics1_v2.txt"};
    //TFile* inFile[2];
    TChain* T[2];
    for (int i=0; i<2; i++){
        T[i] = new TChain("T");
        std::ifstream fileList(inName[i]);
        if (!fileList.is_open()) {
            std::cerr << "Error: Unable to open file list: " << inName[i] << std::endl;
            return;
        }
        std::string line;
        while (std::getline(fileList, line)) {
            if (!line.empty()) {
                T[i]->Add(line.c_str());
            }
        }
        fileList.close();
    }

    TH1D* h_Strip[2][2][2];//[stripmap][layer][u/v]
    TH1D* h_clusterSize[2][2][2];//[stripmap][layer][u/v]
    TH2D* h_2dHitMap[2][2];//[stripmap][layer]
    int nStips[2] = {1600,3200};
    for (int i=0; i<2; i++){
        for (int layer=0; layer<2; layer++){
            for (int uv=0; uv<2; uv++){
                h_Strip[i][layer][uv] = new TH1D(TString::Format("h_Strip_%d_%d_%d",i,layer,uv),TString::Format("h_Strip_%d_%d_%d",i,layer,uv), nStips[uv],0,nStips[uv]);
                h_clusterSize[i][layer][uv] = new TH1D(TString::Format("h_clusterSize_%d_%d_%d",i,layer,uv),TString::Format("h_clusterSize_%d_%d_%d",i,layer,uv), 10,0.5,10.5);
                T[i]->Draw(TString::Format("H.gem.m%i.strip.istrip>>h_Strip_%d_%d_%d",layer, i, layer, uv), TString::Format("H.gem.m%i.strip.IsV==%i",layer, uv), "goff");
                T[i]->Draw(TString::Format("H.gem.clust.nstrip>>h_clusterSize_%d_%d_%d", i, layer, uv), TString::Format("H.gem.clust.layer==%d&&H.gem.clust.axis==%d",layer, uv), "goff");
            }
            h_2dHitMap[i][layer] = new TH2D(TString::Format("h_2dHitMap_%d_%d",i,layer), TString::Format("h_2dHitMap_%d_%d",i,layer), nStips[0],20, 120,nStips[1],-35,35);
            T[i]->Draw(TString::Format("H.gem.sp.posY:H.gem.sp.posX>>h_2dHitMap_%d_%d", i, layer), TString::Format("H.gem.sp.layer==%i",layer), "goff");
        }
    }
    TFile* outFile = new TFile("cmp_stripMap_H.root","RECREATE");
    outFile->cd();
    TCanvas* c1 = new TCanvas("c1","c1",1200,800);
    c1->Divide(2,2);
    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    c2->Divide(2,2);
    for (int i=0; i<2; i++){
        for (int layer=0; layer<2; layer++){
            for (int uv=0; uv<2; uv++){
                c1->cd(layer*2+uv+1);
                h_Strip[i][layer][uv]->SetLineColor(i==0?kRed:kBlue);
                h_Strip[i][layer][uv]->GetXaxis()->SetTitle("Strip Number");
                h_Strip[i][layer][uv]->GetYaxis()->SetTitle("Counts");
                h_Strip[i][layer][uv]->SetTitle(TString::Format("Layer %d %s strips", layer+1, uv==0?"U":"V"));
                if (i==0){
                    h_Strip[i][layer][uv]->DrawClone();
                }else{
                    h_Strip[i][layer][uv]->DrawClone("SAME");
                }
                h_Strip[i][layer][uv]->Write();
                c2->cd(layer*2+uv+1);
                h_clusterSize[i][layer][uv]->SetLineColor(i==0?kRed:kBlue);
                h_clusterSize[i][layer][uv]->GetXaxis()->SetTitle("Cluster Size (number of strips)");
                h_clusterSize[i][layer][uv]->GetYaxis()->SetTitle("Counts");
                h_clusterSize[i][layer][uv]->SetTitle(TString::Format("Layer %d %s cluster size", layer+1, uv==0?"U":"V"));
                if (i==0){
                    h_clusterSize[i][layer][uv]->DrawClone();
                }else{
                    h_clusterSize[i][layer][uv]->DrawClone("SAME");
                }
                h_clusterSize[i][layer][uv]->Write();
            }
            h_2dHitMap[i][layer]->SetTitle(TString::Format("2D Hit Map Layer %d", layer+1));
            h_2dHitMap[i][layer]->GetXaxis()->SetTitle("sp.posX");
            h_2dHitMap[i][layer]->GetYaxis()->SetTitle("sp.posY");
            h_2dHitMap[i][layer]->Write();
        }
    }
    outFile->cd();
    c1->Write();
    c2->Write();
    for (int i=0; i<2; i++){
        delete T[i];
    }
    outFile->Close();

}