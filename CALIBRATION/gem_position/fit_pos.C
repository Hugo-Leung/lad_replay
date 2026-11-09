#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include <map>

const int MAX_DATA  = 10000;
const int maxTracks = 500;
const int nPlanes   = 5;

const double d0_cut               = 20.0;
const double delta_trans_cut[2]   = {-100, 100.0};
const double delta_long_cut[2]    = {-50.0, 50.0};
const double plane_theta[nPlanes] = {150.0, 150.0, 127.0, 127.0, 104.0}; // Angle in degrees
const double plane_r[nPlanes]     = {615.0, 655.6, 523.0, 563.6, 615.0}; // Radius of the second point

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
  double dtheta[2];  // GEM dPhi offsets
};
struct EventData{
    TVector3 gem_hit_pos_local[2];
    TVector3 hodoHit;
    TVector3 targetPos;
};

TVector3 GetHodoHitPosition(const int paddle, const int plane) {

  double dTrans = (5 - paddle) * 22.0;
  // Define the plane using three points
  TVector3 p_hit(plane_r[plane] * cos((plane_theta[plane] - 90) * TMath::DegToRad()), 0,
                 -plane_r[plane] * sin((plane_theta[plane] - 90) * TMath::DegToRad()));
  p_hit = p_hit + TVector3(-dTrans * cos((180 - plane_theta[plane]) * TMath::DegToRad()), 0,
                           -dTrans * sin((180 - plane_theta[plane]) * TMath::DegToRad()));

  return p_hit;
}

TVector3 GetGEMHitPosition(const double x_loc, const double y_loc, const double dx, const double dy, const double r,
                           const double gem_theta, const double gem_dTheta) {
  // Define the plane using three points
  TVector3 p_hit(r * cos((gem_theta - 90) * TMath::DegToRad()), 0, -r * sin((gem_theta - 90) * TMath::DegToRad()));
  p_hit = p_hit + TVector3(-(x_loc + dx) * cos((180 - gem_theta-gem_dTheta) * TMath::DegToRad()), dy,
                           -(x_loc + dx) * sin((180 - gem_theta-gem_dTheta) * TMath::DegToRad()));

  return p_hit;
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
    double resolution[]={20, 0.01, 0.01};//assume 20 cm resolution for hodos, 0.1 mm for GEM
    
    for (int i = 0; i < points.size(); i++) {
        TVector3 dist = diffFromTrack(pTarget, direction, points[i]);
        //chisq += dist * dist; // Sum of squared distances
        //ignore y for now
        chisq += (dist.X()*dist.X() + dist.Z()*dist.Z())/resolution[i]*resolution[i];
    }
    return chisq;
}

TVector3 FitPoints(TVector3* p, TVector3 pTarget, int nPoints, double &fittedChisq) {
    //fit a line through the points, and return the direction vector, must originate from pTarget
    //minimise the chisq function using Tminuit
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimum->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
    minimum->SetTolerance(1e-6);
    minimum->SetPrintLevel(0);
    //double initialParams[3] = {1.0, 0.0, 0.0}; // Initial guess for direction vector (unit vector along x-axis)
    double initialParams[] = {127*TMath::DegToRad()};// Initial guess for direction vector relative to z axis
    ROOT::Math::Functor f([=](const double* params) {
        TVector3 dir(1.0*TMath::Sin(params[0]),0, 1.0*TMath::Cos(params[0]));
        dir = dir.Unit(); // Ensure it's a unit vector
        return chisq_perTrack(pTarget, dir, std::vector<TVector3>(p, p + nPoints));
    }, 1);
    minimum->SetFunction(f);
    minimum->SetVariable(0, "dir_x", initialParams[0], 0.1);
    minimum->Minimize();
    const double* fittedParams = minimum->X();
    fittedChisq = minimum->MinValue();
    //
    TVector3 fittedDirection(1.0*TMath::Sin(fittedParams[0]), 0, 1.0*TMath::Cos(fittedParams[0]));
    fittedDirection = fittedDirection.Unit(); // Ensure it's a unit vector
    return fittedDirection;
}


//GEMParams init_Guess = {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, -1.10238}, {0.0, 0.0},{-8.01005,-8.78698}};
double gemParams_init[6][2] = {{127.018, 127.0}, {0.0, 0.0}, {78.503, 96.6958}, {0.00135288, -1.10759}, {0.0, 0.0},{-8.35324,-8.78311}};
//double gemParams_init[6][2] = {{127.596, 127}, {0.0, 0.0}, {77.8078, 97.0138}, {0.158371, 0.674957}, {0.0, 0.0},{-6,-6}};
double gemParams_tmp[6][2] = {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, -1.10238}, {0.0, 0.0},{-8.01005,-8.78698}};



double total_chisq(GEMParams gemPos, std::vector<EventData> eventDataList){
    double cost = 0.0;
    int nPoints = 0;
    for(const auto& track : eventDataList){
        //Get GEM hit positions with current parameters
        TVector3 p_gem_hit0 = GetGEMHitPosition(track.gem_hit_pos_local[0].X(), track.gem_hit_pos_local[0].Y(),
                                             gemPos.dx[0], gemPos.dy[0], gemPos.r[0], gemPos.theta[0], gemPos.dtheta[0]);
        TVector3 p_gem_hit1 = GetGEMHitPosition(track.gem_hit_pos_local[1].X(), track.gem_hit_pos_local[1].Y(),
                                             gemPos.dx[1], gemPos.dy[1], gemPos.r[1], gemPos.theta[1], gemPos.dtheta[1]);
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
            TVector3 tmpDir = FitPoints(hitPoints, p_react, 3, tmpChisq);
            if(tmpChisq<fittedChisq){
                fittedChisq=tmpChisq;
                fittedDir=tmpDir;
                iVertex=i;
            }
        }
        delete[] hitPoints;
        if (fittedChisq>50){// ignore these tracks for now, they are likely to be randoms
            continue;
        }
        cost += fittedChisq;
        nPoints++;
        if (nPoints%1000==0){
            std::cout << nPoints << " fitted \r" << std::flush;
        }
    }
    return cost/nPoints; 

}

GEMParams toParam(double array[6][2]){
    GEMParams p;
    for(int i =0; i<2; i++){
        p.theta[i]=array[0][i];
        p.phi[i]=array[1][i];
        p.r[i]=array[2][i];
        p.dx[i]=array[3][i];
        p.dy[i]=array[4][i];
        p.dtheta[i]=array[5][i];
    }
    return p;
}

void fit_pos(){
    TString inList = "filelist_optics1.txt";
    int run_number = 0; // Run number to process
    TString spec = "H";// Spectrometer name (H or P)

    std::ifstream fileList(inList);
    if (!fileList.is_open()) {
        std::cerr << "Error: Unable to open file list: " << inList << std::endl;
        return;
    }
    TChain chain("T");
    std::string line;
    while (std::getline(fileList, line)) {
        if (!line.empty()) {
            chain.Add(line.c_str());
        }
    }
    TChain* T = &chain;
    fileList.close();
    TString outputFileName = "gem_pos_fitting_optics.root";


    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
        return; 
    }
    //Define variables
    Double_t trk_d0[MAX_DATA];
    Double_t trk_d0_good[MAX_DATA];
    Double_t trk_projz[MAX_DATA];
    Double_t trk_x[2][MAX_DATA], trk_y[2][MAX_DATA], trk_z[2][MAX_DATA];
    Double_t trk_x_local[2][MAX_DATA], trk_y_local[2][MAX_DATA];
    Double_t kin_trackID_0[MAX_DATA], kin_trackID_1[MAX_DATA];
    Double_t kin_plane_0[MAX_DATA], kin_plane_1[MAX_DATA];
    Double_t kin_paddle_0[MAX_DATA], kin_paddle_1[MAX_DATA];
    Double_t kin_hittime_0[MAX_DATA], kin_hittime_1[MAX_DATA];
    Double_t kin_hittheta_0[MAX_DATA], kin_hittheta_1[MAX_DATA];
    Double_t kin_hitphi_0[MAX_DATA], kin_hitphi_1[MAX_DATA];
    Double_t kin_hitedep_0[MAX_DATA], kin_hitedep_1[MAX_DATA];
    Double_t kin_deltapostrans_0[MAX_DATA], kin_deltapostrans_1[MAX_DATA];
    Double_t kin_deltaposlong_0[MAX_DATA], kin_deltaposlong_1[MAX_DATA];
    Double_t kin_dp, kin_npeSum;
    Int_t nTracks, nGoodHits, nGoodHits1;
    Double_t vertex_x, vertex_y, vertex_z;

    T->SetBranchAddress("Ndata."+spec+".gem.trk.d0", &nTracks);
    T->SetBranchAddress("Ndata."+spec+".ladkin.goodhit_trackid_0", &nGoodHits);
    T->SetBranchAddress("Ndata."+spec+".ladkin.goodhit_trackid_1", &nGoodHits1);
    T->SetBranchAddress(spec+".gem.trk.d0", &trk_d0);
    T->SetBranchAddress(spec+".gem.trk.d0_good", &trk_d0_good);
    T->SetBranchAddress(spec+".gem.trk.projz", &trk_projz);
    T->SetBranchAddress(spec+".gem.trk.x1", &trk_x[0]);
    T->SetBranchAddress(spec+".gem.trk.y1", &trk_y[0]);
    T->SetBranchAddress(spec+".gem.trk.z1", &trk_z[0]);
    T->SetBranchAddress(spec+".gem.trk.x2", &trk_x[1]);
    T->SetBranchAddress(spec+".gem.trk.y2", &trk_y[1]);
    T->SetBranchAddress(spec+".gem.trk.z2", &trk_z[1]);
    T->SetBranchAddress(spec+".gem.trk.x1_local", &trk_x_local[0]);
    T->SetBranchAddress(spec+".gem.trk.y1_local", &trk_y_local[0]);
    T->SetBranchAddress(spec+".gem.trk.x2_local", &trk_x_local[1]);
    T->SetBranchAddress(spec+".gem.trk.y2_local", &trk_y_local[1]);
    T->SetBranchAddress(spec+".ladkin.goodhit_trackid_0", &kin_trackID_0);
    T->SetBranchAddress(spec+".ladkin.goodhit_trackid_1", &kin_trackID_1);
    T->SetBranchAddress(spec+".ladkin.goodhit_plane_0", &kin_plane_0);
    T->SetBranchAddress(spec+".ladkin.goodhit_plane_1", &kin_plane_1);
    T->SetBranchAddress(spec+".ladkin.goodhit_paddle_0", &kin_paddle_0);
    T->SetBranchAddress(spec+".ladkin.goodhit_paddle_1", &kin_paddle_1);
    T->SetBranchAddress(spec+".ladkin.goodhit_hittime_0", &kin_hittime_0);
    T->SetBranchAddress(spec+".ladkin.goodhit_hittime_1", &kin_hittime_1);
    T->SetBranchAddress(spec+".ladkin.goodhit_hittheta_0", &kin_hittheta_0);
    T->SetBranchAddress(spec+".ladkin.goodhit_hittheta_1", &kin_hittheta_1);
    T->SetBranchAddress(spec+".ladkin.goodhit_hitphi_0", &kin_hitphi_0);
    T->SetBranchAddress(spec+".ladkin.goodhit_hitphi_1", &kin_hitphi_1);
    T->SetBranchAddress(spec+".ladkin.goodhit_hitedep_0", &kin_hitedep_0);
    T->SetBranchAddress(spec+".ladkin.goodhit_hitedep_1", &kin_hitedep_1);
    //T->SetBranchAddress(spec+".ladkin.goodhit_deltapostrans_0", &kin_deltapostrans_0);
    //T->SetBranchAddress(spec+".ladkin.goodhit_deltapostrans_1", &kin_deltapostrans_1);
    //T->SetBranchAddress(spec+".ladkin.goodhit_deltaposlong_0", &kin_deltaposlong_0);
    //T->SetBranchAddress(spec+".ladkin.goodhit_deltaposlong_1", &kin_deltaposlong_1);
    T->SetBranchAddress(spec+".react.x", &vertex_x);
    T->SetBranchAddress(spec+".react.y", &vertex_y);
    T->SetBranchAddress(spec+".react.z", &vertex_z);
    T->SetBranchAddress(spec+".gtr.dp", &kin_dp);
    T->SetBranchAddress(spec+".cer.npeSum", &kin_npeSum);
    //T->SetBranchAddress(spec+".ngcer.npeSum", &kin_npeSum);// for SHMS

    int nEntries = T->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;
    // Loop through the tree entries and create a vector of EventData
    std::vector<EventData> eventDataList;
    for (Long64_t i = 0; i < nEntries; ++i) {
        T->GetEntry(i);
        // Skip processing if there are too many tracks
        if (nTracks > maxTracks) {
            continue;
        }
        
        if(kin_npeSum<2 || kin_dp<-10.0 || kin_dp>10.0 || vertex_z>20 || vertex_z<-20) continue; // basic cuts for HMS

        // Adjust vertex_z to the nearest fixed z position for optics 1 target
        for (int i=0; i<nFixedz; i++){
            if(vertex_z>=z_range[i] && vertex_z<z_range[i+1]){
                vertex_z = fixedz[i];
                break;
            }
        }
         for (int j = 0; j < nGoodHits; ++j){
            int trackID0 = static_cast<int>(kin_trackID_0[j]);
            int trackID1 = static_cast<int>(kin_trackID_1[j]);
            if ((trackID0 < 0 || trackID0 >= nTracks) && (trackID1 < 0 || trackID1 >= nTracks)) {
                continue;
            }
            int planeID, paddleID;
            int trackID=-1;
            TVector3 p_hit;
            if (trackID0 >= 0 && trackID0 < nTracks) {
                planeID = static_cast<int>(kin_plane_0[j]);
                paddleID = static_cast<int>(kin_paddle_0[j]);
                p_hit = GetHodoHitPosition(paddleID, planeID);
                trackID=trackID0;
            }else if (trackID1 >= 0 && trackID1 < nTracks) {
                planeID = static_cast<int>(kin_plane_1[j]);
                paddleID = static_cast<int>(kin_paddle_1[j]);
                p_hit = GetHodoHitPosition(paddleID, planeID);
                trackID=trackID1;
            }
            TVector3 p_react(vertex_x, vertex_y, vertex_z);
            EventData eventData;
            eventData.hodoHit = p_hit;
            eventData.targetPos = p_react;
            // Get GEM hit positions in local coordinates
            for (int gemIndex = 0; gemIndex < 2; ++gemIndex) {
                double x_loc = (gemIndex == 0) ? trk_x_local[0][trackID] : trk_x_local[1][trackID];
                double y_loc = (gemIndex == 0) ? trk_y_local[0][trackID] : trk_y_local[1][trackID];
                TVector3 gem_hit_pos (x_loc, y_loc, 0.0); // Z=0 in local coordinates
                eventData.gem_hit_pos_local[gemIndex] = gem_hit_pos;
            }
            eventDataList.push_back(eventData);
         }
         if(i%1000==0){
            std::cout<<"\r filling "<<i<<"/"<<nEntries<<std::flush;
         }
    }
    std::cout<< "finished filling hit list"<<std::endl;
    // Now perform fitting to optimize GEM positions
    // Initial guess for GEM parameters
    GEMParams gemParams = toParam(gemParams_init);

    outputFile->cd();
    
    for (int itr=0; itr<12*5; itr++){
        int index=itr%12;
        if (index/2==1 || index/2==4){
            std::cout<<"skipping "<<itr<<std::endl;
            continue; //skiping phi and y for now
        }
        std::cout<<"itr: "<<itr<<std::endl;
        TGraph* gr_chisq = new TGraph();
        TGraph* gr_dTheta =  new TGraph();
    
        ROOT::Math::Functor f([=](const double* params) {
            for (int j=0; j<6;j++){
                for (int i=0; i<2; i++){
                    gemParams_tmp[j][i]=gemParams_init[j][i];
                }
            }
            int index=(int)params[1];
            gemParams_tmp[index/2][index%2]=params[0];
            GEMParams gemParams = toParam(gemParams_tmp);
            double chisq = total_chisq(gemParams, eventDataList);
            gr_chisq->SetPoint(gr_chisq->GetN(),gr_chisq->GetN(), chisq);
            gr_dTheta->SetPoint(gr_dTheta->GetN(), gr_dTheta->GetN(), params[0]);
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
        const double minChisq = minimum->MinValue();
        std::cout << "Min chisq: " << minChisq <<" at: " << fittedParams[0] << std::endl;
        gemParams_init[index/2][index%2]=fittedParams[0];

        gr_chisq->SetMarkerStyle(21);
        gr_chisq->SetName(TString::Format("chisq_itr_%i",itr));
        gr_chisq->Write();
        gr_dTheta->SetMarkerStyle(22);
        gr_dTheta->SetName(TString::Format("Param_%i_%i_itr_%i",index/2,index%2,itr));
        gr_dTheta->Write();
    }
    //clean up
    outputFile->Write();
    outputFile->Close();
    for (int i=0; i<6; i++){
        std::cout<<"{";
        for (int j=0; j<2; j++){
            std::cout<<gemParams_init[i][j]<<" ";
        }
        std::cout<<"}"<<endl;

    }
    std::cout << "Fitting completed. Results saved to " << outputFileName << std::endl;
}