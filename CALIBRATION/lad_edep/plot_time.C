void plot_time(){
    const int N_PLANES                 = 5;
    const int N_PADDLES                = 11;
    TString inName = "edep_time.dat";
    TString outName = "edep_time.root";
    TFile *file = TFile::Open(outName, "RECREATE");
    if (!file || file->IsZombie()) {
        cout << "Error: Cannot create output ROOT file." << endl;
        return;
    }
    TTree *T = new TTree("T", "Edep Time Tree");
    T->ReadFile(inName);
    TGraphErrors *g[2];
    int plane, paddle;
    double time, time_err;
    const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
    T->SetBranchAddress("plane", &plane);
    T->SetBranchAddress("bar", &paddle);
    T->SetBranchAddress("t", &time);
    T->SetBranchAddress("et", &time_err);
    int nEntries = T->GetEntries();
    for (int i=0; i<2; ++i) {
        g[i] = new TGraphErrors();
        g[i]->SetName(Form("g_time_plane_%s", plane_names[2*i+1].c_str()));
        g[i]->SetTitle(Form("Time vs Paddle for Plane %s", plane_names[2*i+1].c_str()));
        g[i]->GetXaxis()->SetTitle("Paddle");
        g[i]->GetYaxis()->SetTitle("Time (ns)");
    }
    for (int i=0; i<nEntries; ++i) {
        T->GetEntry(i);
        if (plane % 2 == 0) continue; // Only process odd planes
        int idx = (plane - 1) / 2; // Convert plane to index
        g[idx]->SetPoint(g[idx]->GetN(), paddle, time);
        g[idx]->SetPointError(g[idx]->GetN()-1, 0, time_err);
    }
    TCanvas *c = new TCanvas("c_time", "Time vs Paddle", 800, 600);
    c->Divide(1, 2);
    for (int i=0; i<2; ++i) {
        c->cd(i+1);
        g[i]->SetMarkerStyle(20);
        g[i]->SetMarkerColor(i+1);
        g[i]->Draw("AP");
        g[i]->Fit("pol1", "Q");
        g[i]->GetFunction("pol1")->SetLineColor(i+1);
        g[i]->GetFunction("pol1")->SetLineWidth(2);
        g[i]->GetFunction("pol1")->Draw("same");
    }
    c->Update();
    c->Write();
    delete c;
    delete T;
    for (int i=0; i<2; ++i) {
        delete g[i];
    }
    file->Close();
    delete file;
    cout << "Time plots saved to " << outName << endl;
    return;
}