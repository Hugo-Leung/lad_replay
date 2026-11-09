void test(){
    TChain* T = new TChain("T");
    std::string line="/volatile/hallc/c-lad/cleung/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_23463_0_0_1_80000.root";
    T->Add(line.c_str());

    double triggerBits=0;

    T->SetBranchAddress("g.trigbits", &triggerBits);
    int nEntries = T->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        T->GetEntry(entry);
        if(entry>=1e2) break;
        UInt_t tBits = static_cast<UInt_t>(triggerBits);
        std::cout<<entry<<" "<<tBits<<" "<<(tBits&8)<<" "<<(tBits&2)<<std::endl;
    }
}