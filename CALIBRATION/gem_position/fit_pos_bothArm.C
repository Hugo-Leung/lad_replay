#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include <map>

const int MAX_DATA  = 10000;
const int maxTracks = 10000;
const int nPlanes   = 5;

const double d0_cut               = 20.0;
const double delta_trans_cut[2]   = {-100, 100.0};
const double delta_long_cut[2]    = {-50.0, 50.0};
const double plane_theta[nPlanes] = {149.6825077, 149.8727149, 126.7569325, 126.844879, 103.9683928}; // Angle in degrees
const double plane_r[nPlanes]     = {618.8125006, 658.5516251, 526.1277344, 568.079082, 614.5200935}; // Radius of the second point
const double plane_y[nPlanes]    = {3.429, 1.463, 2.232, 1.538, 1.538}; // y position of the plane

const int nFixedz = 3;
// const double fixedz[nFixedz] = {0.0}; // Fixed z positions for the planes
const double fixedz[nFixedz] = {-10.0, 0.0, 10.0}; // Fixed z positions for the planes
// const double z_range[nFixedz + 1] = {-15.0, -5.0, 5.0, 15.0}; // Tolerance range for the fixed z positions (low,
// high)
const double z_range[nFixedz + 1] = {-20.0, -5, 5.0, 20.0}; // Tolerance range for the fixed z positions (low,
// high)

struct GEMParams {
  double theta[2]; // GEM angles in degrees
  double phi[2];   // GEM phi angles in degrees
  double r[2];     // GEM radii
  double dx[2];    // GEM dx offsets
  double dy[2];    // GEM dy offsets
  double dphi[2];  // GEM dPhi offsets
  double dtheta[2];  // GEM dPhi offsets
  double dRoll[2];  // GEM dRoll offsets
};
struct EventData{
    TVector3 gem_hit_pos_local[2];
    TVector3 hodoHit;
    TVector3 targetPos;
    double dir[2]={127.*TMath::DegToRad(),0};
};

TVector3 GetHodoHitPosition(const int paddle, const int plane) {

  double dTrans = (5 - paddle) * 22.0;
  // Define the plane using three points
  TVector3 p_hit(plane_r[plane] * cos((plane_theta[plane] - 90) * TMath::DegToRad()), plane_y[plane],
                 -plane_r[plane] * sin((plane_theta[plane] - 90) * TMath::DegToRad()));
  p_hit = p_hit + TVector3(-dTrans * cos((180 - plane_theta[plane]) * TMath::DegToRad()), 0,
                           -dTrans * sin((180 - plane_theta[plane]) * TMath::DegToRad()));

  return p_hit;
}

TVector3 GetGEMHitPosition(const double x_loc, const double y_loc, const double dx, const double dy, const double r,
                           const double gem_theta, const double gem_dTheta, const double gem_phi, const double gem_dPhi, const double gem_dRoll) {
    TVector3 p_loc(x_loc, y_loc, 0.0); // Local coordinates
    p_loc.RotateX(gem_dPhi * TMath::DegToRad()); // Rotate by dPhi around X-axis
    p_loc.RotateY(gem_dTheta * TMath::DegToRad()); // Rotate by dTheta around Y-axis
    p_loc.RotateZ(gem_dRoll * TMath::DegToRad()); // Rotate by dRoll around Z-axis
    p_loc.SetX(p_loc.X() + dx);
    p_loc.SetY(p_loc.Y() + dy);
    p_loc.SetZ(p_loc.Z()+r);
    p_loc.RotateY(gem_theta * TMath::DegToRad()); // Rotate by theta around Y-axis
    p_loc.RotateZ(gem_phi * TMath::DegToRad()); // Rotate by theta around Y-axis
    return p_loc;
}

TVector3 diffFromTrack (TVector3 pTarget, TVector3 direction, TVector3 point) {
    // Calculate the vector from the target position to the point
    TVector3 toPoint = point - pTarget;

    // Project the toPoint vector onto the direction vector
    double projectionLength = toPoint.Dot(direction)/direction.Mag2();//direction should be a unit vector, but just in case
    TVector3 projection     = projectionLength * direction;

    // Calculate the perpendicular distance from the point to the line
    TVector3 perpendicular = toPoint - projection;
    //double distance       = perpendicular.Mag();

    return perpendicular;
}

double chisq_perTrack(TVector3 pTarget, TVector3 direction, std::vector<TVector3> points) {
    double chisq = 0.0;
    double resolution[]={10, 0.1, 0.1};//assume 10 cm resolution for hodos, 0.1 mm for GEM
    double yhodo_res = 10; // 2.5 cm resolution in y for hodo
    for (int i = 0; i < points.size(); i++) {
        TVector3 dist = diffFromTrack(pTarget, direction, points[i]);
        //chisq += dist * dist; // Sum of squared distances
        if (i==0){
            if ((dist.X()*dist.X() + dist.Z()*dist.Z())<resolution[i]*resolution[i]) // if hodo residual is less than resolution, ignore
                chisq += 0;
            else{
                chisq += (dist.X()*dist.X() + dist.Z()*dist.Z())/(resolution[i]*resolution[i]);
            }
            chisq += (dist.Y()*dist.Y())/(yhodo_res*yhodo_res);
        }
        else
        chisq += (dist.X()*dist.X() +dist.Y()*dist.Y()+ dist.Z()*dist.Z())/(resolution[i]*resolution[i]);
    }
    return chisq;
}

TVector3 FitPoints(TVector3* p, TVector3 pTarget, int nPoints, double &fittedChisq, double initialParams[2]=nullptr) {
    //fit a line through the points, and return the direction vector, must originate from pTarget
    //minimise the chisq function using Tminuit
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimum->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
    minimum->SetTolerance(1e-6);
    minimum->SetPrintLevel(0);
    //double initialParams[3] = {1.0, 0.0, 0.0}; // Initial guess for direction vector (unit vector along x-axis)
    if (initialParams==nullptr)
        initialParams = new double[2]{127*TMath::DegToRad(),0};// theta and phi for direction vector relative to z axis
    ROOT::Math::Functor f([=](const double* params) {
        TVector3 dir(1.0*TMath::Sin(params[0])*TMath::Cos(params[1]), 1.0*TMath::Sin(params[0])*TMath::Sin(params[1]), 1.0*TMath::Cos(params[0]));
        dir = dir.Unit(); // Ensure it's a unit vector
        return chisq_perTrack(pTarget, dir, std::vector<TVector3>(p, p + nPoints));
    }, 2);
    minimum->SetFunction(f);
    minimum->SetVariable(0, "theta", initialParams[0], 0.1);
    minimum->SetVariable(1, "phi", initialParams[1], 0.1);
    minimum->Minimize();
    const double* fittedParams = minimum->X();
    fittedChisq = minimum->MinValue();
    //
    TVector3 fittedDirection(1.0*TMath::Sin(fittedParams[0])*TMath::Cos(fittedParams[1]), 1.0*TMath::Sin(fittedParams[0])*TMath::Sin(fittedParams[1]), 1.0*TMath::Cos(fittedParams[0]));
    fittedDirection = fittedDirection.Unit(); // Ensure it's a unit vector
    delete minimum;
    return fittedDirection;
}


//double gemParams_init[8][2] = {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 0.0}, {0.0, 0.0},{0,0},{0,0}, {0,0}};
double gemParams_init[8][2] = {{127.008, 127 }, {0.057755 ,0.000285756}, {85.9126,  105.913}, { 0.0390553 ,-0.0337049}, {1.33479, 1.33553}, {-1.1876,  -0.322357 }, { 1.26238,  0.0888885}, {-0.049502, -0.0762915}}; //16.2635
//double gemParams_init[8][2] = {{126.693, 126.969}, {1.01428, 0.816717}, {81.5718, 99.873}, {0.945316, 1.51709}, {0.319978, 0.205513}, {-13.1795, -15.9727}, {2.00119, 0.521598}, {0.301969, 0.333304}};
double gemParams_tmp[8][2] = {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 0.0}, {0.0, 0.0},{0,0},{-8.01005,-8.78698}, {0,0}};



double total_chisq(GEMParams gemPos, std::vector<EventData> eventDataList){
    double cost = 0.0;
    int nPoints = 0;
    for( auto& track : eventDataList){
        //Get GEM hit positions with current parameters
        TVector3 p_gem_hit0 = GetGEMHitPosition(track.gem_hit_pos_local[0].X(), track.gem_hit_pos_local[0].Y(),
                                             gemPos.dx[0], gemPos.dy[0], gemPos.r[0], gemPos.theta[0], gemPos.dtheta[0], gemPos.phi[0], gemPos.dphi[0], gemPos.dRoll[0]);
        TVector3 p_gem_hit1 = GetGEMHitPosition(track.gem_hit_pos_local[1].X(), track.gem_hit_pos_local[1].Y(),
                                             gemPos.dx[1], gemPos.dy[1], gemPos.r[1], gemPos.theta[1], gemPos.dtheta[1], gemPos.phi[1], gemPos.dphi[1], gemPos.dRoll[1]);
        TVector3 p_hit = track.hodoHit;
        TVector3 p_react = track.targetPos;
        TVector3* hitPoints = new TVector3[3];
        hitPoints[0] = p_hit;
        hitPoints[1] = p_gem_hit0;
        hitPoints[2] = p_gem_hit1;
        double fittedChisq = 1e10;
        TVector3 fittedDir;
        int iVertex=-1;
        for (int i=0; i<nFixedz; i++){
            double tmpChisq=0.0;
            p_react.SetZ(fixedz[i]);
            TVector3 tmpDir = FitPoints(hitPoints, p_react, 3, tmpChisq, track.dir);
            if(tmpChisq<fittedChisq){
                fittedChisq=tmpChisq;
                fittedDir=tmpDir;
                iVertex=i;
            }
        }
        delete[] hitPoints;
        track.dir[0]=TMath::ACos(fittedDir.Z());
        track.dir[1]=TMath::ATan2(fittedDir.Y(), fittedDir.X());
        if (fittedChisq>10){// ignore these tracks for now, they are likely to be randoms
            continue;
        }
        cost += fittedChisq;
        nPoints++;
        if (nPoints%1000==0){
            std::cout << nPoints << " fitted \r" << std::flush;
        }
    }
    std::cout << nPoints << " fitted; chi2/nTrack: "<<cost/nPoints<<"\r" << std::flush;
    return cost/nPoints; 

}

GEMParams toParam(double array[8][2]){
    GEMParams p;
    for(int i =0; i<2; i++){
        p.theta[i]=array[0][i];
        p.phi[i]=array[1][i];
        p.r[i]=array[2][i];
        p.dx[i]=array[3][i];
        p.dy[i]=array[4][i];
        p.dphi[i]=array[5][i];
        p.dtheta[i]=array[6][i];
        p.dRoll[i]=array[7][i];
    }
    return p;
}

void fit_pos_bothArm(int iParam_start){
    TString inList = "filelist_optics1_17deg.txt";
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
    TString outputFileName = "gem_pos_fitting_optics1_17deg_bothArms.root";

    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
        return; 
    }

    //Define variables for HMS and SHMS
    Double_t trk_d0[2][MAX_DATA];
    Double_t trk_d0_good[2][MAX_DATA];
    Double_t trk_projz[2][MAX_DATA];
    Double_t trk_x[2][2][MAX_DATA], trk_y[2][2][MAX_DATA], trk_z[2][2][MAX_DATA];
    Double_t trk_x_local[2][2][MAX_DATA], trk_y_local[2][2][MAX_DATA];
    Double_t kin_trackID_0[2][MAX_DATA], kin_trackID_1[2][MAX_DATA];
    Double_t kin_plane_0[2][MAX_DATA], kin_plane_1[2][MAX_DATA];
    Double_t kin_paddle_0[2][MAX_DATA], kin_paddle_1[2][MAX_DATA];
    Double_t kin_hittime_0[2][MAX_DATA], kin_hittime_1[2][MAX_DATA];
    Double_t kin_hittheta_0[2][MAX_DATA], kin_hittheta_1[2][MAX_DATA];
    Double_t kin_hitphi_0[2][MAX_DATA], kin_hitphi_1[2][MAX_DATA];
    Double_t kin_hitedep_0[2][MAX_DATA], kin_hitedep_1[2][MAX_DATA];
    Double_t kin_deltapostrans_0[2][MAX_DATA], kin_deltapostrans_1[2][MAX_DATA];
    Double_t kin_deltaposlong_0[2][MAX_DATA], kin_deltaposlong_1[2][MAX_DATA];
    Double_t kin_dp[2], kin_npeSum[2];
    Int_t nTracks[2], nGoodHits[2], nGoodHits1[2];
    Double_t kin_yhodo_0[2][MAX_DATA], kin_yhodo_1[2][MAX_DATA];
    Double_t vertex_x[2], vertex_y[2], vertex_z[2];



    UInt_t triggerBits=0;
    Double_t triggerBits_float=0;
    UInt_t triggerMask[]={8,2};// HMS Ereal, SHMS Ereal
    
    T->SetBranchAddress("g.trigbits", &triggerBits_float);
    //set branch addresses for HMS and SHMS
    for (int iSpec=0; iSpec<nSpecs; iSpec++){
        T->SetBranchAddress("Ndata."+spec[iSpec]+".gem.trk.d0", &nTracks[iSpec]);
        T->SetBranchAddress("Ndata."+spec[iSpec]+".ladkin.goodhit_trackid_0", &nGoodHits[iSpec]);
        T->SetBranchAddress("Ndata."+spec[iSpec]+".ladkin.goodhit_trackid_1", &nGoodHits1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.d0", &trk_d0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.d0_good", &trk_d0_good[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.projz", &trk_projz[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.x1", &trk_x[iSpec][0]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.y1", &trk_y[iSpec][0]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.z1", &trk_z[iSpec][0]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.x2", &trk_x[iSpec][1]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.y2", &trk_y[iSpec][1]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.z2", &trk_z[iSpec][1]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.x1_local", &trk_x_local[iSpec][0]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.y1_local", &trk_y_local[iSpec][0]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.x2_local", &trk_x_local[iSpec][1]);
        T->SetBranchAddress(spec[iSpec]+".gem.trk.y2_local", &trk_y_local[iSpec][1]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_trackid_0", &kin_trackID_0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_trackid_1", &kin_trackID_1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_plane_0", &kin_plane_0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_plane_1", &kin_plane_1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_paddle_0", &kin_paddle_0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_paddle_1", &kin_paddle_1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_ypos_0", &kin_yhodo_0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_ypos_1", &kin_yhodo_1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hittime_0", &kin_hittime_0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hittime_1", &kin_hittime_1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hittheta_0", &kin_hittheta_0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hittheta_1", &kin_hittheta_1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hitphi_0", &kin_hitphi_0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hitphi_1", &kin_hitphi_1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hitedep_0", &kin_hitedep_0[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".ladkin.goodhit_hitedep_1", &kin_hitedep_1[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".react.x", &vertex_x[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".react.y", &vertex_y[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".react.z", &vertex_z[iSpec]);
        T->SetBranchAddress(spec[iSpec]+".gtr.dp", &kin_dp[iSpec]);
        if (iSpec==0)
            T->SetBranchAddress(spec[iSpec]+".cer.npeSum", &kin_npeSum[iSpec]);
        else
            T->SetBranchAddress(spec[iSpec]+".ngcer.npeSum", &kin_npeSum[iSpec]);
    }

    int nEntries = T->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;

    //Collect event data
    std::vector<EventData> eventDataList;
    for (int entry = 0; entry < nEntries; ++entry) {
        T->GetEntry(entry);
        // Skip processing if there are too many tracks
        int iSpec=-1;
        triggerBits = static_cast<UInt_t>(triggerBits_float);
        //Check trigger
        if ((triggerBits&triggerMask[0])!=0){//HMS
            iSpec=0;
            //if(kin_npeSum[iSpec]<2 || kin_dp[iSpec]<-10.0 || kin_dp[iSpec]>10.0 || vertex_z[iSpec]>20 || vertex_z[iSpec]<-20) continue; // basic cuts for HMS
        }
        if ((triggerBits&triggerMask[1])!=0){//SHMS
            iSpec=1;
            //if(kin_npeSum[iSpec]<6 || kin_dp[iSpec]<-15.0 || kin_dp[iSpec]>24.0 || vertex_z[iSpec]>20 || vertex_z[iSpec]<-20) continue; // basic cuts for SHMS
        }
        if (iSpec==-1){
            std::cout<<"no valid spectrometer triggered"<<std::endl;
            std::cout<<entry<<" "<<triggerBits<<" "<<(triggerBits&triggerMask[0])<<" "<<(triggerBits&triggerMask[1])<<std::endl;
            return ;
            //continue; // no valid spectrometer triggered

        }

        if (nTracks[iSpec] > maxTracks) {
            continue;
        }
        EventData event;

        TVector3 p_react(vertex_x[iSpec], vertex_y[iSpec], vertex_z[iSpec]);

        // Adjust vertex_z to the nearest fixed z position for optics 1 target
        for (int i=0; i<nFixedz; i++){
            if(p_react.Z()>=z_range[i] && p_react.Z()<z_range[i+1]){
                p_react.SetZ( fixedz[i]);
                break;
            }
        }

         for (int j = 0; j < nGoodHits[iSpec]; ++j){
            int trackID0 = static_cast<int>(kin_trackID_0[iSpec][j]);
            int trackID1 = static_cast<int>(kin_trackID_1[iSpec][j]);
            if ((trackID0 < 0 || trackID0 >= nTracks[iSpec]) && (trackID1 < 0 || trackID1 >= nTracks[iSpec])) {
                continue;
            }
            
            int planeID, paddleID;
            int trackID=-1;
            TVector3 p_hit;
            if (trackID0 >= 0 && trackID0 < nTracks[iSpec]) {
                planeID = static_cast<int>(kin_plane_0[iSpec][j]);
                paddleID = static_cast<int>(kin_paddle_0[iSpec][j]);
                p_hit = GetHodoHitPosition(paddleID, planeID);
                if (kin_yhodo_0[iSpec][j]<-1000 || kin_yhodo_0[iSpec][j]>1000){//invalid ypos from LADKIN, skip
                    continue;
                }
                p_hit.SetY(kin_yhodo_0[iSpec][j]); // use ypos from LADKIN
                trackID=trackID0;
            }else if (trackID1 >= 0 && trackID1 < nTracks[iSpec]) {
                planeID = static_cast<int>(kin_plane_1[iSpec][j]);
                paddleID = static_cast<int>(kin_paddle_1[iSpec][j]);
                p_hit = GetHodoHitPosition(paddleID, planeID);
                if (kin_yhodo_1[iSpec][j]<-1000 || kin_yhodo_1[iSpec][j]>1000){//invalid ypos from LADKIN, skip
                    continue;
                }
                p_hit.SetY(kin_yhodo_1[iSpec][j]); // use ypos from LADKIN
                trackID=trackID1;
            }
            EventData eventData;
            eventData.hodoHit = p_hit;
            eventData.targetPos = p_react;
            // Get GEM hit positions in local coordinates
            for (int gemIndex = 0; gemIndex < 2; ++gemIndex) {
                double x_loc = (gemIndex == 0) ? trk_x_local[iSpec][0][trackID] : trk_x_local[iSpec][1][trackID];
                double y_loc = (gemIndex == 0) ? trk_y_local[iSpec][0][trackID] : trk_y_local[iSpec][1][trackID];
                TVector3 gem_hit_pos (x_loc, y_loc, 0.0); // Z=0 in local coordinates
                eventData.gem_hit_pos_local[gemIndex] = gem_hit_pos;
            }
            eventDataList.push_back(eventData);
         }
         if(entry%1000==0){
            std::cout<<"\r filling "<<entry<<"/"<<nEntries<<"; size(eventDataList): "<<eventDataList.size()<<std::flush;
         }

    }
    std::cout<< "finished filling hit list"<<std::endl;
    // Now perform fitting to optimize GEM positions
    // Initial guess for GEM parameters
    delete T;
    outputFile->cd();
    GEMParams gemParams = toParam(gemParams_init);
    TGraph* gr_chisq = new TGraph();
    TGraph* gr_param[16];
    for (int i=0; i<16; i++){
        gr_param[i] = new TGraph();
        gr_param[i]->SetName(TString::Format("param_%d", i).Data());
    }

    
    for (int itr=0; itr<16*10; itr++){
        int index=(itr+iParam_start)%16;
        std::cout<<"itr: "<<itr<<std::endl;
        /*if (index/2<5||index/2>8){
            std::cout<<"skip "<<index <<" for now"<<std::endl;
            gr_chisq->SetPoint(itr, itr, total_chisq(toParam(gemParams_init), eventDataList));
            gr_param[index]->SetPoint(itr, itr, gemParams_init[index/2][index%2]);
            continue;
        }*/
    
        ROOT::Math::Functor f([=](const double* params) {
            for (int j=0; j<8;j++){
                for (int i=0; i<2; i++){
                    gemParams_tmp[j][i]=gemParams_init[j][i];
                }
            }
            int index=(int)params[1];
            gemParams_tmp[index/2][index%2]=params[0];
            GEMParams gemParams = toParam(gemParams_tmp);
            double chisq = total_chisq(gemParams, eventDataList);
            return chisq;
        }, 1);
        ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimum->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
        minimum->SetTolerance(1e-4);
        minimum->SetFunction(f);
        minimum->SetVariable(0, TString::Format("param_%i_%i ",index/2,index%2).Data(),gemParams_init[index/2][index%2] , 0.1);
        minimum->SetVariable(1,"index",index,0);
        minimum->FixVariable(1);
        minimum->Minimize();

        const double* fittedParams = minimum->X();
        const double* errorEst = minimum->Errors();
        const double minChisq = minimum->MinValue();
        std::cout << "Min chisq: " << minChisq <<" at: " << fittedParams[0] <<" +/- "<<errorEst[0]<< std::endl;
        gemParams_init[index/2][index%2]=fittedParams[0];
        gr_chisq->SetPoint(itr, itr, minChisq);
        gr_param[index]->SetPoint(itr, itr, fittedParams[0]);
        delete minimum;
        
    }
    //clean up
    for (int i=0; i<16; i++){
        gr_param[i]->Write();
    }
    gr_chisq->Write();
    outputFile->Write();
    outputFile->Close();
    for (int i=0; i<8; i++){
        std::cout<<"{";
        for (int j=0; j<2; j++){
            std::cout<<gemParams_init[i][j]<<" ";
        }
        std::cout<<"}"<<endl;

    }
    std::cout << "Fitting completed. Results saved to " << outputFileName << std::endl;
}
