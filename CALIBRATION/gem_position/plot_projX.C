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

std::map<int, GEMParams> runToGEMParams = {
    {0, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{0,0}} },
    {1, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{-6.30242,-6.30242}} }
};

const double z_targ     = 0;   // z position of the target plane in cm

TVector3 gemPos[2];
TVector3 gemPlaneDir[2];


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

TVector3 GEMWindow(const int paddle, const int plane, const int gemIndex, const TVector3 vertex){
    TVector3 hodoHit = GetHodoHitPosition(paddle, plane);
    TVector3 dir = (hodoHit - vertex);
    //calculate intersection of this line with the GEM plane
    //plane point
    TVector3 p0 = gemPos[gemIndex];
    double t = dir.Dot(p0 - vertex) / (p0-vertex).Dot(p0-vertex);
    TVector3 p_hit = vertex + t*dir;
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
    double resolution[]={10, 0.01, 0.01};//assume 10 cm resolution for hodos, 0.1 mm for GEM
    
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
    std::cout << "Fitted chisq: " << fittedChisq << std::endl;
    //
    TVector3 fittedDirection(1.0*TMath::Sin(fittedParams[0]), 0, 1.0*TMath::Cos(fittedParams[0]));
    fittedDirection = fittedDirection.Unit(); // Ensure it's a unit vector
    return fittedDirection;
}



void plot_projX(){
    int runNumber = 0; // Run number to process

    TString inList = "filelist_optics1.txt";
    int run_number = runNumber; // Run number to process
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
    TString outputFileName = TString::Format("gem_pos_fixing_%i.root",runNumber);

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
    // Define histograms
    TH1D* h_gem_hit_along_plane[2];
    TH1D* h_vertex_z = new TH1D("h_vertex_z", "Vertex Z Distribution;Vertex Z (cm);Counts", 100, -20, 20);
    TH1D* h_nGoodTracks = new TH1D("h_nGoodTracks", "Number of Good Tracks;Number of Good Tracks;Counts", 100, 0, maxTracks+10);
    TH1D* h_nGoodHits = new TH1D("h_nGoodHits", "Number of Good Hits;Number of Good Hits;Counts", 20, 0, 20);
    TH1D* h_hodo_hit_along_plane[2];
    TH2D* h_gem_diff_from_track[2];
    TH2D* h_hodo_diff_from_track[2];
    TH2D* h_gem_hit_local[2];
    TH2D* h_hodo_hit_position = new TH2D("h_hodo_hit_position", "Hodo hit position;X (cm);z (cm)", 70, 0, 700, 70, -700, 0);
    TH1D* h_hodo_diff_target;
    h_hodo_diff_target = new TH1D("h_hodo_diff_target", "Hodo hit diff from target;Diff from target (cm);Counts", 100, -50, 50);
    TH1D* h_planeID= new TH1D("h_planeID", "Hodo plane ID;Plane ID;Counts", 6, -0.5, 5.5);
    TH1D* h_fitted_chisq = new TH1D("h_fitted_chisq", "Fitted Chisq;Chisq;Counts", 100, 0, 500);
    TH1D* h_fitted_dir = new TH1D("h_fitted_dir", "Fitted Direction Angle;Angle (rad);Counts", 100, 0, 2*TMath::Pi());
    for(int i=0; i<2; i++){
        h_gem_hit_local[i] = new TH2D(Form("h_gem_hit_local_gem%d", i+1), Form("GEM %d local hit position;Local X (cm);Local Y (cm)", i+1), 130, -65, 65, 120, -60, 60);
        h_gem_hit_along_plane[i] = new TH1D(Form("h_gem_hit_along_plane_gem%d", i+1), Form("GEM %d hit along plane;Position along plane (cm);Counts", i+1), 130, -65, 65);
        h_hodo_hit_along_plane[i] = new TH1D(Form("h_hodo_hit_along_plane_hodo%d", i+1), Form("Hodo hit along plane for Hodo Plane %d;Position along plane (cm);Counts", i+1), 100, -200, 200);
        h_gem_diff_from_track[i] = new TH2D(Form("h_gem_diff_from_track_gem%d", i+1), Form("GEM %d diff from track vs position along plane;Position along plane (cm);Diff from track (cm)", i+1), 130, -65, 65, 100, -30, 30);
        h_hodo_diff_from_track[i] = new TH2D(Form("h_hodo_diff_from_track_hodo%d", i+1), Form("Hodo Plane %d diff from track vs position along plane;Position along plane (cm);Diff from track (cm)", i+1), 100, -200, 200, 100, -30, 30);
    }

    gemPos[0] = TVector3(runToGEMParams[run_number].r[0] * cos((runToGEMParams[run_number].theta[0] - 90) * TMath::DegToRad()), 0,
                         -runToGEMParams[run_number].r[0] * sin((runToGEMParams[run_number].theta[0] - 90) * TMath::DegToRad()));
    gemPos[1] = TVector3(runToGEMParams[run_number].r[1] * cos((runToGEMParams[run_number].theta[1] - 90) * TMath::DegToRad()), 0,
                         -runToGEMParams[run_number].r[1] * sin((runToGEMParams[run_number].theta[1] - 90) * TMath::DegToRad()));

    gemPlaneDir[0] = TVector3(-cos((180 - runToGEMParams[run_number].theta[0]-runToGEMParams[run_number].dtheta[0]) * TMath::DegToRad()), 0,
                             -sin((180 - runToGEMParams[run_number].theta[0]-runToGEMParams[run_number].dtheta[0]) * TMath::DegToRad()));
    gemPlaneDir[1] = TVector3(-cos((180 - runToGEMParams[run_number].theta[1]-runToGEMParams[run_number].dtheta[1]) * TMath::DegToRad()), 0,
                             -sin((180 - runToGEMParams[run_number].theta[1]-runToGEMParams[run_number].dtheta[1]) * TMath::DegToRad()));

    for(int i =0; i <nEntries; i++){
        T->GetEntry(i);
        //std::cout << "vetex_z: "<<vertex_z<<"; nTracks: " << nTracks << std::endl;
        if (nTracks > maxTracks) {
            continue;
        }
        if(i%1000 == 0){
            std::cout << "Processing entry " << i << " / " << nEntries << "\n" << std::flush;
        }
        //if(i>1e5)   break; // for testing
        //if(kin_npeSum<6 || kin_dp<-15.0 || kin_dp>24.0 || vertex_z>10 || vertex_z<-10) continue; // basic cuts for SHMS
        if(kin_npeSum<2 || kin_dp<-10.0 || kin_dp>10.0 || vertex_z>20 || vertex_z<-20) continue; // basic cuts for HMS
        h_vertex_z->Fill(vertex_z);
        //vertex_z = 0;// For calibration data, set vertex_z to 0, single foil target

        for (int i=0; i<nFixedz; i++){
            if(vertex_z>=z_range[i] && vertex_z<z_range[i+1]){
                vertex_z = fixedz[i];
                break;
            }
        }

        h_nGoodTracks->Fill(nTracks);
        h_nGoodHits->Fill(nGoodHits);
        // Create a vector to label tracks as good or not
        std::vector<bool> goodTrack(nTracks, false);
        // Create vectors to store hodoscope paddle, plane, and y position
        std::vector<int> hodoPaddle(nTracks, -1);
        std::vector<int> hodoPlane(nTracks, -1);
        std::vector<double> hodoYPos(nTracks, 0.0);
        // Loop through all hodo hits
        for (int j = 0; j < nGoodHits; ++j) {
            // Get the track ID associated with the hodo hit
            int trackID0 = static_cast<int>(kin_trackID_0[j]);
            int trackID1 = static_cast<int>(kin_trackID_1[j]);
            /*if(trackID != trackID1){
                std::cout << "Warning: Mismatched track IDs for front and back hodoscope hits! front: "<<trackID<<"; back: "<<trackID1<<std::endl;
            } */

            // Ensure the track ID is within bounds
            if ((trackID0 < 0 || trackID0 >= nTracks) && (trackID1 < 0 || trackID1 >= nTracks)) {
                continue;
            }
            h_planeID->Fill(kin_plane_0[j]);
            h_planeID->Fill(kin_plane_1[j]);
            // get hodo hit position
            int planeID, paddleID;
            int trackID=-1;
            TVector3 p_hit;
            if (trackID0 >= 0 && trackID0 < nTracks) {
                trackID=trackID0;
                planeID = static_cast<int>(kin_plane_0[j]);
                paddleID = static_cast<int>(kin_paddle_0[j]);
                p_hit = GetHodoHitPosition(paddleID, planeID);
            }else if (trackID1 >= 0 && trackID1 < nTracks) {
                trackID=trackID1;
                planeID = static_cast<int>(kin_plane_1[j]);
                paddleID = static_cast<int>(kin_paddle_1[j]);
                p_hit = GetHodoHitPosition(paddleID, planeID);
            }
            
            TVector3 p_react(vertex_x, vertex_y, vertex_z);
            TVector3 gemHits[2];
            for(int gemIndex=0; gemIndex<2; gemIndex++){
                if (trk_x_local[gemIndex][trackID]==trk_y_local[gemIndex][trackID]&&trk_y_local[gemIndex][trackID]==0){
                    std::cout<<trackID<<endl;
                    return;

                }
                gemHits[gemIndex] = GetGEMHitPosition(trk_x_local[gemIndex][trackID], trk_y_local[gemIndex][trackID],
                                                     runToGEMParams[run_number].dx[gemIndex],
                                                     runToGEMParams[run_number].dy[gemIndex],
                                                     runToGEMParams[run_number].r[gemIndex],
                                                     runToGEMParams[run_number].theta[gemIndex],runToGEMParams[run_number].dtheta[gemIndex]);
                h_gem_hit_local[gemIndex]->Fill(trk_x_local[gemIndex][trackID], trk_y_local[gemIndex][trackID]);
            }
            TVector3* hitPoints = new TVector3[3];
            hitPoints[0] = p_hit;
            hitPoints[1] = gemHits[0];
            hitPoints[2] = gemHits[1];
            double fittedChisq = 0.0;
            TVector3 fittedDir = FitPoints(hitPoints, p_react, 3, fittedChisq);
            if (fittedChisq > 100) {
                delete[] hitPoints;
                continue;
            }
            h_fitted_chisq->Fill(fittedChisq);
            h_fitted_dir->Fill(atan2(fittedDir.X(), fittedDir.Z()));
            h_hodo_hit_position->Fill(p_hit.X(), p_hit.Z());
            delete[] hitPoints;

            h_gem_diff_from_track[0]->Fill(trk_x_local[0][trackID], diffFromTrack(p_react, fittedDir, gemHits[0]).Dot(gemPlaneDir[0]));
            h_gem_diff_from_track[1]->Fill(trk_x_local[1][trackID], diffFromTrack(p_react, fittedDir, gemHits[1]).Dot(gemPlaneDir[1]));

            
        }
    }
    // End Event Loop
    std::cout << "\nProcessing completed." << std::endl;
    // Save histograms to the output file
    /*outputFile->cd();
    h_vertex_z->Write();
    h_hodo_diff_target->Write();
    h_hodo_hit_position->Write();
    h_nGoodTracks->Write();
    h_nGoodHits->Write();
    for(int i=0; i<2; i++){
        h_gem_hit_along_plane[i]->Write();
        h_hodo_hit_along_plane[i]->Write();
        h_gem_diff_from_track[i]->Write();
        h_hodo_diff_from_track[i]->Write();
    }*/
    // write and close the file
    outputFile->Write();
    outputFile->Close();
    std::cout << "Output written to " << outputFileName << std::endl;
    return;

}
