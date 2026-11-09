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

TVector3 projectToTarget(TVector3 pTarget, TVector3 gemHit1, TVector3 gemHit2){
    TVector3 dir = gemHit2-gemHit1;
    double t = (pTarget.X()-gemHit1.X())/dir.X();
    TVector3 projection = gemHit1+t*dir;
    return projection;

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


std::map<int, GEMParams> runToGEMParams = {
    {0, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{0,0}} },
    {1, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{-6.30242,-6.30242}} },
    {2, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{-12.9884,-12.9884}} },
    {3, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{-8.01005,-8.01005}} },
    {4, {{127.902, 127}, {0.0, 0.0}, {77.571, 97.0138}, {0.158371, 1.11008}, {0.0, 0.0},{0,0}} },
 //   {5, {{127.596, 127}, {0.0, 0.0}, {77.8078, 97.0138}, {0.158371, 0.674957}, {0.0, 0.0},{1.96605,-2.09132}} }
    {5, {{127.902, 127}, {0.0, 0.0}, {77.571, 97.0138}, {0.158371, 1.11008}, {0.0, 0.0},{-8.00945, -8.78839}} },
    {6, {{127.019, 127}, {0.0, 0.0}, {78.525, 96.7154}, {0.00143581, -1.10771}, {0.0, 0.0}, {-8.35201, -8.76677}} }

};

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
        double fittedChisq = 0.0;
        TVector3 fittedDir = FitPoints(hitPoints, p_react, 3, fittedChisq);
        delete[] hitPoints;
        if (fittedChisq>100){// ignore these tracks for now, they are likely to be randoms
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


void plot_residual(int iParam){
    TString inList = "filelist_optics1.txt";
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
    TString outputFileName = TString::Format("gem_pos_fixed_optics_%i.root",iParam);


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

    int nEntries = T->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;
    // Loop through the tree entries and create a vector of EventData
    std::vector<EventData> eventDataList;

    TH2D* h_gem_hit_2D[2];
    TH1D* h_vertex_z = new TH1D("h_vertex_z", "Vertex Z Distribution;Vertex Z (cm);Counts", 100, -20, 20);
    TH1D* h_nGoodTracks = new TH1D("h_nGoodTracks", "Number of Good Tracks;Number of Good Tracks;Counts", 100, 0, maxTracks+10);
    TH1D* h_nGoodHits = new TH1D("h_nGoodHits", "Number of Good Hits;Number of Good Hits;Counts", 20, 0, 20);
    TH2D* h_gem_diff_from_track[2];
    TH2D* h_gem_hit_local[2];
    for(int i=0; i<2; i++){
        h_gem_hit_local[i] = new TH2D(Form("h_gem_hit_local_gem%d", i+1), Form("GEM %d local hit position;Local X (cm);Local Y (cm)", i+1), 130, -65, 65, 120, -60, 60);
        h_gem_hit_2D[i] = new TH2D(Form("h_gem_hit_2d_gem%d", i+1), Form("GEM %d hit Position; X (cm); z (cm)", i+1), 120, 0, 120, 140, -120, 20);
        h_gem_diff_from_track[i] = new TH2D(Form("h_gem_diff_from_track_gem%d", i+1), Form("GEM %d diff from track vs position along plane;Position along plane (cm);Diff from track (cm)", i+1), 130, -65, 65, 100, -30, 30);
    }
    TH2D* h_hodo_hit_position = new TH2D("h_hodo_hit_position", "Hodo hit position;X (cm);z (cm)", 70, 0, 700, 70, -700, 0);
    TH1D* h_fitted_chisq = new TH1D("h_fitted_chisq", "Fitted Chisq;Chisq;Counts", 100, 0, 100);
    TH1D* h_fitted_dir = new TH1D("h_fitted_dir", "Fitted Direction Angle;Angle (rad);Counts", 100, 1, 3);

    TH2D* h_vertex_gem = new TH2D("h_vertex_gem", "Vertex Z; gem Vertex Z(cm);electron vertex Z(cm)",100,-20,20, 100,-20,20);
    
    GEMParams gemPos = runToGEMParams[iParam];
    TVector3 gemPlaneDir[2];
    gemPlaneDir[0] = TVector3(-cos((180 - gemPos.theta[0]-gemPos.dtheta[0]) * TMath::DegToRad()), 0,
                             -sin((180 - gemPos.theta[0]-gemPos.dtheta[0]) * TMath::DegToRad()));
    gemPlaneDir[1] = TVector3(-cos((180 - gemPos.theta[1]-gemPos.dtheta[1]) * TMath::DegToRad()), 0,
                             -sin((180 - gemPos.theta[1]-gemPos.dtheta[1]) * TMath::DegToRad()));
    for (Long64_t i = 0; i < nEntries; ++i) {
        T->GetEntry(i);
        // Skip processing if there are too many tracks
        if (nTracks > maxTracks) {
            continue;
        }
        if(kin_npeSum<2 || kin_dp<-10.0 || kin_dp>10.0 || vertex_z>20 || vertex_z<-20) continue; // basic cuts for HMS

        h_vertex_z->Fill(vertex_z);
        double eVertex_z=vertex_z;
        // Adjust vertex_z to the nearest fixed z position for optics 1 target
        for (int i=0; i<nFixedz; i++){
            if(vertex_z>=z_range[i] && vertex_z<z_range[i+1]){
                vertex_z = fixedz[i];
                break;
            }
        }
        h_nGoodTracks->Fill(nTracks);
        h_nGoodHits->Fill(nGoodHits);
         for (int j = 0; j < nGoodHits; ++j){
            int trackID0 = static_cast<int>(kin_trackID_0[j]);
            int trackID1 = static_cast<int>(kin_trackID_1[j]);
            if ((trackID0 < 0 || trackID0 >= nTracks) && (trackID1 < 0 || trackID1 >= nTracks)) {
                continue;
            }
            int planeID, paddleID;
            TVector3 p_hit;
            int trackID=-1;
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
            h_hodo_hit_position->Fill(p_hit.X(), p_hit.Z());
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
                h_gem_hit_local[gemIndex]->Fill(x_loc, y_loc);
            }
            eventDataList.push_back(eventData);


            TVector3 p_gem_hit0 = GetGEMHitPosition(eventData.gem_hit_pos_local[0].X(), eventData.gem_hit_pos_local[0].Y(),
                                                 gemPos.dx[0], gemPos.dy[0], gemPos.r[0], gemPos.theta[0], gemPos.dtheta[0]);
            TVector3 p_gem_hit1 = GetGEMHitPosition(eventData.gem_hit_pos_local[1].X(), eventData.gem_hit_pos_local[1].Y(),
                                                 gemPos.dx[1], gemPos.dy[1], gemPos.r[1], gemPos.theta[1], gemPos.dtheta[1]);
            TVector3* hitPoints = new TVector3[3];
            hitPoints[0] = p_hit;
            hitPoints[1] = p_gem_hit0;
            hitPoints[2] = p_gem_hit1;

            double fittedChisq = 1e10;
            TVector3 fittedDir;
            int iVertex=-1;
            for (int j=0; j<nFixedz; j++){
                double tmpChisq=0.0;
                p_react.SetZ(fixedz[j]);
                TVector3 tmpDir = FitPoints(hitPoints, p_react, 3, tmpChisq);
                if(tmpChisq<fittedChisq){
                    fittedChisq=tmpChisq;
                    fittedDir=tmpDir;
                    iVertex=j;
                }
            }
            delete[] hitPoints;
            if (fittedChisq>50){// ignore these tracks for now, they are likely to be randoms
                continue;
            }
            p_react.SetZ(fixedz[iVertex]);
            h_fitted_chisq->Fill(fittedChisq);
            h_fitted_dir->Fill(atan2(fittedDir.X(), fittedDir.Z()));
            h_hodo_hit_position->Fill(p_hit.X(), p_hit.Z());
            h_gem_hit_2D[0]->Fill(p_gem_hit0.X(),p_gem_hit0.Z());
            h_gem_hit_2D[1]->Fill(p_gem_hit1.X(),p_gem_hit1.Z());
            h_gem_diff_from_track[0]->Fill(eventData.gem_hit_pos_local[0].X(), diffFromTrack(p_react, fittedDir, p_gem_hit0).Dot(gemPlaneDir[0]));
            h_gem_diff_from_track[1]->Fill(eventData.gem_hit_pos_local[1].X(), diffFromTrack(p_react, fittedDir, p_gem_hit1).Dot(gemPlaneDir[1]));
            h_vertex_gem->Fill(projectToTarget(p_react,p_gem_hit0,p_gem_hit1).Z(),eVertex_z);
         }
        if (i%1000==0){
            std::cout<<"\rprocessiing "<<i<<"/"<<nEntries<<std::flush;
        }
    }

    TH1D* h_vertex_proj= h_vertex_gem->ProjectionX();
    TF1* fit = new TF1("fit_3guas","gaus(0)+gaus(3)+gaus(6)+[9]",-20,20);
    Double_t fit_p[]={20, -10,2, 20, 0, 2, 14, 10, 2,7};
    fit->SetParameters(fit_p);
    fit->SetParNames("Norm1","Mean1", "SD1", "Norm2","Mean2", "SD2", "Norm3","Mean3", "SD3","offset");
    fit->SetParLimits(0,0,1e5);
    fit->SetParLimits(3,0,1e5);
    fit->SetParLimits(6,0,1e5);
    fit->SetParLimits(1,-20,0);
    fit->SetParLimits(4,-10,10);
    fit->SetParLimits(7,0,20);
    fit->SetParLimits(2,0.1,10);
    fit->SetParLimits(5,0.1,10);
    fit->SetParLimits(8,0.1,10);    
    TFitResultPtr fitResult =h_vertex_proj->Fit(fit,"S");
    
    TF1* g1 = new TF1("g1", "gaus(0)",-20,20);
    g1->SetParameter(0,fitResult->Parameter(0));
    g1->SetParameter(1,fitResult->Parameter(1));
    g1->SetParameter(2,fitResult->Parameter(2));

    g1->SetLineColor(kBlue);
    h_vertex_proj->GetListOfFunctions()->Add(g1);

    TF1* g2 = new TF1("g2", "gaus(0)",-20,20);
    g2->SetParameter(0,fitResult->Parameter(3));
    g2->SetParameter(1,fitResult->Parameter(4));
    g2->SetParameter(2,fitResult->Parameter(5));

    g2->SetLineColor(kGreen+2);
    h_vertex_proj->GetListOfFunctions()->Add(g2);
    
    TF1* g3 = new TF1("g1", "gaus(0)",-20,20);
    g3->SetParameter(0,fitResult->Parameter(6));
    g3->SetParameter(1,fitResult->Parameter(7));
    g3->SetParameter(2,fitResult->Parameter(8));

    g3->SetLineColor(kBlack);
    h_vertex_proj->GetListOfFunctions()->Add(g3);
    
    TF1* c = new TF1("const", "[0]",-20,20);
    c->SetParameter(0,fitResult->Parameter(9));

    c->SetLineColor(kCyan);
    h_vertex_proj->GetListOfFunctions()->Add(c);

    double maxY = h_vertex_proj->GetMaximum();

    for (int i=0; i<nFixedz; i++){
        TLine* l = new TLine(fixedz[i],0,fixedz[i],maxY);
        l->SetLineColor(14);
        h_vertex_proj->GetListOfFunctions()->Add(l); 
        TLine* l1 = new TLine (fixedz[i], h_vertex_gem->GetYaxis()->GetXmin(), fixedz[i], h_vertex_gem->GetYaxis()->GetXmax());
        l1->SetLineColor(14);
        h_vertex_gem->GetListOfFunctions()->Add(l1);
        for (int j=0; j<nFixedz; j++){
            TLine* l2 = new TLine (h_vertex_gem->GetXaxis()->GetXmin(), fixedz[j], h_vertex_gem->GetXaxis()->GetXmax(), fixedz[j]);
            l2->SetLineColor(14);
            h_vertex_gem->GetListOfFunctions()->Add(l2);
        }

    }
    h_vertex_z->Fit(fit);
    

    //clean up
    outputFile->Write();
    outputFile->Close();
    std::cout << "Proccssing completed. Results saved to " << outputFileName << std::endl;
}