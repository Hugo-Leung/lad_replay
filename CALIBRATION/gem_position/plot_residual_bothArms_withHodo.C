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
    p_loc.RotateZ(gem_phi * TMath::DegToRad()); // Rotate by phi around Z-axis
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
    double yhodo_res = 10; // 1.5 cm resolution in y for hodo
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

TVector3 FitPoints(TVector3* p, TVector3 pTarget, int nPoints, double &fittedChisq, double &fittedZ) {
    //fit a line through the points, and return the direction vector, must originate from pTarget
    //minimise the chisq function using Tminuit
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimum->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
    minimum->SetTolerance(1e-6);
    minimum->SetPrintLevel(0);
    //double initialParams[3] = {1.0, 0.0, 0.0}; // Initial guess for direction vector (unit vector along x-axis)
    double initialParams[3] = {127*TMath::DegToRad(),0,0};// theta and phi for direction vector relative to z axis and z vertex
    ROOT::Math::Functor f([=](const double* params) {
        TVector3 pTarget_mod = pTarget;
        pTarget_mod.SetZ(params[2]);
        TVector3 dir(1.0*TMath::Sin(params[0])*TMath::Cos(params[1]), 1.0*TMath::Sin(params[0])*TMath::Sin(params[1]), 1.0*TMath::Cos(params[0]));
        dir = dir.Unit(); // Ensure it's a unit vector
        return chisq_perTrack(pTarget_mod, dir, std::vector<TVector3>(p, p + nPoints));
    }, 3);
    minimum->SetFunction(f);
    minimum->SetVariable(0, "theta", initialParams[0], 0.1);
    minimum->SetVariable(1, "phi", initialParams[1], 0.1);
    minimum->SetVariable(2, "z_vertex", initialParams[2], 0.1);
    minimum->Minimize();
    const double* fittedParams = minimum->X();
    fittedChisq = minimum->MinValue();
    
    TVector3 fittedDirection(1.0*TMath::Sin(fittedParams[0])*TMath::Cos(fittedParams[1]), 1.0*TMath::Sin(fittedParams[0])*TMath::Sin(fittedParams[1]), 1.0*TMath::Cos(fittedParams[0]));
    fittedDirection = fittedDirection.Unit(); // Ensure it's a unit vector
    fittedZ = fittedParams[2];
    delete minimum;
    return fittedDirection;
}

TVector3 projectToTarget(TVector3 pTarget, TVector3 gemHit1, TVector3 gemHit2){
    TVector3 dir = gemHit2-gemHit1;
    dir = dir.Unit();
    // closest approch to pTarget from gemHit1 along dir using only x and y
    double t = ((pTarget.X()-gemHit1.X())*dir.X() + (pTarget.Y()-gemHit1.Y())*dir.Y())/(dir.X()*dir.X() + dir.Y()*dir.Y());
    TVector3 projection = gemHit1+t*dir;
    return projection;

}


double gemParams_init[8][2] = {{127.902, 127}, {0.0, 0.0}, {77.571, 99.3051}, {0.0, 0.0}, {0.0, 0.0},{0., 0.},{-8.00945, -7.35828},{0., 0.}};
std::map<int, GEMParams> runToGEMParams = {
    {0, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{0., 0.},{0,0},{0., 0.}} },
    {1, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{0., 0.},{-0.430832, -1.72714},{0., 0.}} },
    {2, {{127.002, 127.0}, {0.0, 0.0}, {77.5821, 95.571}, {0.920331, 1.14758}, {0.0, 0.0},{0., 0.},{-0.432006,  -1.74051},{0., 0.}} },
    {3, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0},{0., 0.},{-8.01005,-8.01005},{0., 0.}} },
    {4, {{127.902, 127}, {0.0, 0.0}, {77.571, 97.0138}, {0.158371, 1.11008}, {0.0, 0.0},{0., 0.},{0,0},{0., 0.}} },
    {5, {{127.902, 127}, {0.0, 0.0}, {77.571, 97.0138}, {0,0}, {0.0, 0.0},{0., 0.},{-8.00945, -8.78839},{0., 0.}} },
    {6, {{126.999, 126.941}, {0.0, 0.0}, {80.8899, 98.067}, {0., 1.4996}, {0., 0.}, {0., 0.},{3.42198, 0.548143},{0., 0.}} },
    {7, {{126.646, 126.969}, {0.959864, 0.827659}, {82.046, 99.9377}, {0.0960709, 1.49865}, {0.322269,0.00871353}, {-8.17174,  -10.3063},{1.88289, 0.631218},{0., 0.}}},
    {8, {{126.646, 126.969}, {0.959315, 0.824642}, {81.6089, 99.9652}, {0.207119, 1.49862}, {0.319874,0.0633616}, {-10.2222,  -10.7726},{1.97517, 0.695669},{0.0157642, 0.0154243}}},
    {9, {{126.668, 126.969}, {0.957231, 0.820631}, {80.7101, 99.8956}, {0.211644, 1.49859}, {0.319198, 0.069482}, { -10.238, -10.8078},{1.97704, 0.378223},{0.298961,0.0215506}}},
    {10, {{126.693, 126.969}, {1.01459, 0.81709}, {81.5717, 99.873}, {0.945236, 1.51709}, {0.319978, 0.205513}, {-13.1794, -15.9727}, {2.00118, 0.521598}, {0.301969, 0.333304}}},
    {11, {{126.808, 127 }, {-0.00173414, 0.587734}, {79.0234, 96.2848 }, {-0.205624, -0.00165544}, {1.33877, 0.259285}, { -19.4683, -20.8754}, {0.802703, -0.189912}, {0.380487, 0.561991}}},
    {12, {{126.809, 127 }, {-0.0002708, 0.622982}, {79.0598, 96.2861}, {-0.140077, -0.0016581}, {1.32588, 0.333158}, { -19.4686, -20.8754}, {0.864875, -0.00104045}, {0.427752, 0.562662}}},
    {13, {{126.819, 127}, {0.00169352, 0.622956}, {79.0849, 96.3009}, {-0.136301, -0.00229047}, {1.3248, 0.332222}, {-19.474, -20.8774}, {0.87989, -0.00436277}, {0.440224, 0.567971}}},
    {22, {{126.61, 126.544}, {0.478144 ,-0.0137268}, {83.4624, 104.104}, { 6.93714e-07 ,-2.97805e-06}, {0.786051,  1.3282}, {-2.96567,1.3604 }, {  1.47474,  0.860781}, {-0.0332988, -1.39128e-05}}},
    {23, {{126.595, 126.545}, {0.538376 , -0.0137268}, {84.7974, 104.104}, { -0.035012 ,-0.000403268}, {0.736886,  1.32821}, {1.3244,1.35708 }, {  0.991763,  0.860786}, {0.247536, -0.0537964}}},
    {24, {{126.586, 126.544}, {0.538037 , -0.0137376}, {84.7974, 103.909}, { -0.033407 ,-0.00029584}, {0.737075,  1.32821}, {1.32763,1.35708 }, {  0.990042,  0.86083}, {0.247884, -0.0534267}}},
    {26,  {{126.586, 126.545}, {0.539155 , -0.0134036}, {84.9544, 104.098}, { -0.0448764 ,-0.0415198}, {0.736875,  1.32821}, {1.62662,1.36715 }, {  0.984768,  0.860878}, {0.247847, -0.0538493}}}, //16.2635
    
    // 17deg
    {30, {{127, 127 }, { 0,0}, {85.271, 105.913}, { -0.059123 ,-0.0415148}, {0.736602,  1.32821}, {-0.3168,  -0.3168 }, { 0.0882,  0.0882}, {-0.0756, -0.0756}}},
    {31, {{127.01, 127 }, {0.0575266,0}, {85.2727, 105.913}, { -0.0590348 ,-0.0333941}, {1.35364,  1.32896}, {-0.369055,  -0.316795 }, {  1.99698,  0.0888122}, {-0.075583, -0.0756209}}},

};

double total_chisq(GEMParams gemPos, std::vector<EventData> eventDataList){
    double cost = 0.0;
    int nPoints = 0;
    for(const auto& track : eventDataList){
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
            double fittedZ=0.0;
            TVector3 tmpDir = FitPoints(hitPoints, p_react, 3, tmpChisq, fittedZ);
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

void plot_residual_bothArms_withHodo(int iParam){
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
    TString outputFileName = TString::Format("gem_pos_fixed_residual_optics_17deg_bothArms_withHodo%i.root",iParam);

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

    GEMParams gemPos = runToGEMParams[iParam];
    TVector3 gemPlaneDir[2][2];//how to compute this?
    for(int i=0; i<2; i++){
        gemPlaneDir[i][0] = TVector3(1,0,0);
    
        gemPlaneDir[i][0].RotateX(gemPos.dphi[i]*TMath::DegToRad());
        gemPlaneDir[i][0].RotateY(gemPos.dtheta[i]*TMath::DegToRad());
        gemPlaneDir[i][0].RotateZ(gemPos.dRoll[i]*TMath::DegToRad());
        gemPlaneDir[i][0].RotateY(gemPos.theta[i]*TMath::DegToRad());
        gemPlaneDir[i][0].RotateZ(gemPos.phi[i]*TMath::DegToRad());
        //
        gemPlaneDir[i][1] = TVector3(0,1,0);
        gemPlaneDir[i][1].RotateX(gemPos.dphi[i]*TMath::DegToRad());
        gemPlaneDir[i][1].RotateY(gemPos.dtheta[i]*TMath::DegToRad());
        gemPlaneDir[i][1].RotateZ(gemPos.dRoll[i]*TMath::DegToRad());
        gemPlaneDir[i][1].RotateY(gemPos.theta[i]*TMath::DegToRad());
        gemPlaneDir[i][1].RotateZ(gemPos.phi[i]*TMath::DegToRad());
    }
    //initialize all the histograms
    TH2D* h_gem_hit_2D[2];
    TH1D* h_vertex_z = new TH1D("h_vertex_z", "Vertex Z Distribution;Vertex Z (cm);Counts", 100, -20, 20);
    TH1D* h_nGoodTracks = new TH1D("h_nGoodTracks", "Number of Good Tracks;Number of Good Tracks;Counts", 100, 0, maxTracks+10);
    TH1D* h_nGoodHits = new TH1D("h_nGoodHits", "Number of Good Hits;Number of Good Hits;Counts", 20, 0, 20);
    TH2D* h_gem_diff_x[2];
    TH2D* h_gem_diff_y[2];
    TH2D* h_gem_diff_xy[2];
    TH2D* h_gem_diff_yx[2];
    TH2D* h_gem_hit_local[2];
    TH1D* h_gem_hit_y[2];
    for(int i=0; i<2; i++){
        h_gem_hit_local[i] = new TH2D(Form("h_gem_hit_local_gem%d", i+1), Form("GEM %d local hit position;Local X (cm);Local Y (cm)", i+1), 130, -65, 65, 120, -60, 60);
        h_gem_hit_2D[i] = new TH2D(Form("h_gem_hit_2d_gem%d", i+1), Form("GEM %d hit Position; X (cm); z (cm)", i+1), 120, 0, 120, 140, -120, 20);
        h_gem_diff_x[i] = new TH2D(Form("h_gem_diff_x_gem%d", i+1), Form("GEM %d diff from track vs position along plane; local x (cm);Diff x (cm)", i+1), 130, -65, 65, 100, -2, 2);
        h_gem_diff_y[i] = new TH2D(Form("h_gem_diff_y_gem%d", i+1), Form("GEM %d diff from track vs position along plane; local y (cm);Diff y (cm)", i+1), 130, -35, 30, 100, -2, 2);
        h_gem_diff_xy[i] = new TH2D(Form("h_gem_diff_xy_gem%d", i+1), Form("GEM %d diff from track vs position along plane; local x (cm);Diff y (cm)", i+1), 130, -65, 65, 100, -2, 2);
        h_gem_diff_yx[i] = new TH2D(Form("h_gem_diff_yx_gem%d", i+1), Form("GEM %d diff from track vs position along plane; local y (cm);Diff x (cm)", i+1), 130, -35, 30, 100, -2, 2);
        h_gem_hit_y[i] = new TH1D(Form("h_gem_hit_y_gem%d", i+1), Form("GEM %d hit Y position;Y (cm);Counts", i+1), 120, -60, 60);
    }
    TH2D* h_hodo_hit_position = new TH2D("h_hodo_hit_position", "Hodo hit position;X (cm);z (cm)", 70, 0, 700, 70, -700, 0);
    TH1D* h_fitted_chisq = new TH1D("h_fitted_chisq", "Fitted Chisq;Chisq;Counts", 100, 0, 25);
    TH1D* h_fitted_dir = new TH1D("h_fitted_dir", "Fitted Direction Angle;Angle (rad);Counts", 100, 1, 3);
    TH1D* h_fitted_dir_y = new TH1D("h_fitted_dir_y", "Fitted Direction Angle;Angle (rad);Counts", 100, -1, 1);

    TH2D* h_vertex_gem = new TH2D("h_vertex_gem", "Vertex Z; gem Vertex Z(cm);electron vertex Z(cm)",100,-20,20, 100,-20,20);
    TH2D* h_vertex_gem_y = new TH2D("h_vertex_gem_y", "Vertex Y; gem Vertex Y(cm);electron vertex Y(cm)",100,-2,2, 100,-0.2,0.2);
    TH1D* h_vertex_gem_y_1D = new TH1D("h_vertex_gem_y_1D", "Vertex Y; gem Vertex Y(cm)",100,-4,4);
    TH2D* h_vertex_gem_x = new TH2D("h_vertex_gem_x", "Vertex X; gem Vertex X(cm);electron vertex X(cm)",100,-2,2, 100,-0.2,0.2);
    TH2D* h_vertex_gem_yz = new TH2D("h_vertex_gem_yz", "Vertex Y vs Z; gem Vertex Z(cm);gem Vertex Y(cm)",100,-20,20, 100,-0.2,0.2);
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
        //    if(kin_npeSum[iSpec]<2 || kin_dp[iSpec]<-10.0 || kin_dp[iSpec]>10.0 || vertex_z[iSpec]>20 || vertex_z[iSpec]<-20) continue; // basic cuts for HMS
        }
        if ((triggerBits&triggerMask[1])!=0){//SHMS
            iSpec=1;
        //    if(kin_npeSum[iSpec]<6 || kin_dp[iSpec]<-15.0 || kin_dp[iSpec]>24.0 || vertex_z[iSpec]>20 || vertex_z[iSpec]<-20) continue; // basic cuts for SHMS
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

        h_vertex_z->Fill(vertex_z[iSpec]);
        double eVertex_z=vertex_z[iSpec];

        TVector3 p_react(vertex_x[iSpec], vertex_y[iSpec], vertex_z[iSpec]);

        h_nGoodTracks->Fill(nTracks[iSpec]);
        h_nGoodHits->Fill(nGoodHits[iSpec]);

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
                p_hit.SetY(kin_yhodo_0[iSpec][j]); // use ypos from LADKIN
                trackID=trackID0;
                if (kin_yhodo_0[iSpec][j]<-1000 || kin_yhodo_0[iSpec][j]>1000){//invalid ypos from LADKIN, skip
                    continue;
                }
            }else if (trackID1 >= 0 && trackID1 < nTracks[iSpec]) {
                planeID = static_cast<int>(kin_plane_1[iSpec][j]);
                paddleID = static_cast<int>(kin_paddle_1[iSpec][j]);
                p_hit = GetHodoHitPosition(paddleID, planeID);
                p_hit.SetY(kin_yhodo_1[iSpec][j]); // use ypos from LADKIN
                trackID=trackID1;
                if (kin_yhodo_1[iSpec][j]<-1000 || kin_yhodo_1[iSpec][j]>1000){//invalid ypos from LADKIN, skip
                    continue;
                }
            }
            h_hodo_hit_position->Fill(p_hit.X(), p_hit.Z());
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
            TVector3 p_gem_hit0 = GetGEMHitPosition(eventData.gem_hit_pos_local[0].X(), eventData.gem_hit_pos_local[0].Y(),
                                                 gemPos.dx[0], gemPos.dy[0], gemPos.r[0], gemPos.theta[0], gemPos.dtheta[0], gemPos.phi[0], gemPos.dphi[0], gemPos.dRoll[0]);
            TVector3 p_gem_hit1 = GetGEMHitPosition(eventData.gem_hit_pos_local[1].X(), eventData.gem_hit_pos_local[1].Y(),
                                                 gemPos.dx[1], gemPos.dy[1], gemPos.r[1], gemPos.theta[1], gemPos.dtheta[1], gemPos.phi[1], gemPos.dphi[1], gemPos.dRoll[1]);
            TVector3* hitPoints = new TVector3[3];
            hitPoints[0] = p_hit;
            hitPoints[1] = p_gem_hit0;
            hitPoints[2] = p_gem_hit1;

            double fittedChisq = 1e10;
            double bestZ=0.0;
            TVector3 fittedDir = FitPoints(hitPoints, p_react, 3, fittedChisq, bestZ);
            delete[] hitPoints;
            if (fittedChisq>10){// ignore these tracks for now, they are likely to be randoms
                continue;
            }
            p_react.SetZ(bestZ);
            h_fitted_chisq->Fill(fittedChisq);
            h_fitted_dir->Fill(acos(fittedDir.Z()));
            h_fitted_dir_y->Fill(atan2(fittedDir.Y(), fittedDir.X()));
            h_hodo_hit_position->Fill(p_hit.X(), p_hit.Z());
            h_gem_hit_2D[0]->Fill(p_gem_hit0.X(),p_gem_hit0.Z());
            h_gem_hit_2D[1]->Fill(p_gem_hit1.X(),p_gem_hit1.Z());
            TVector3 diff0 = diffFromTrack(p_react, fittedDir, p_gem_hit0);
            TVector3 diff1 = diffFromTrack(p_react, fittedDir, p_gem_hit1);
            h_gem_diff_y[0]->Fill(eventData.gem_hit_pos_local[0].Y(), diff0.Dot(gemPlaneDir[0][1]));
            h_gem_diff_x[0]->Fill(eventData.gem_hit_pos_local[0].X(), diff0.Dot(gemPlaneDir[0][0]));
            h_gem_diff_xy[0]->Fill(eventData.gem_hit_pos_local[0].X(), diff0.Dot(gemPlaneDir[0][1]));
            h_gem_diff_yx[0]->Fill(eventData.gem_hit_pos_local[0].Y(), diff0.Dot(gemPlaneDir[0][0]));

            h_gem_diff_y[1]->Fill(eventData.gem_hit_pos_local[1].Y(), diff1.Dot(gemPlaneDir[1][1]));
            h_gem_diff_x[1]->Fill(eventData.gem_hit_pos_local[1].X(), diff1.Dot(gemPlaneDir[1][0]));
            h_gem_diff_xy[1]->Fill(eventData.gem_hit_pos_local[1].X(), diff1.Dot(gemPlaneDir[1][1]));
            h_gem_diff_yx[1]->Fill(eventData.gem_hit_pos_local[1].Y(), diff1.Dot(gemPlaneDir[1][0]));
            h_vertex_gem->Fill(bestZ,eVertex_z);
            h_gem_hit_local[0]->Fill(eventData.gem_hit_pos_local[0].X(), eventData.gem_hit_pos_local[0].Y());
            h_gem_hit_local[1]->Fill(eventData.gem_hit_pos_local[1].X(), eventData.gem_hit_pos_local[1].Y());
            h_vertex_gem_yz->Fill(bestZ, p_react.Y());


         }
         if(entry%1000==0){
            int nTracks = h_fitted_chisq->GetEntries();
            std::cout<<"\r filling "<<entry<<"/"<<nEntries<<"; nTracks: "<<nTracks<<" (estimate: "<<(int)((double)nTracks/entry*nEntries)<<")"<<std::flush;
            /*if (nTracks>300){
                std::cout<<std::endl;
                break;
            }*/
        }

    }
    std::cout<< "finished filling hit list"<<std::endl;
    // Now perform fitting to optimize GEM positions
    // Initial guess for GEM parameters
    delete T;

    outputFile->cd();

    TH1D* h_vertex_proj= h_vertex_gem->ProjectionX();
    TF1* fit = new TF1("fit_3guas","gaus(0)+gaus(3)+gaus(6)",-20,20);
    Double_t fit_p[]={20, -10,2, 20, 0, 2, 14, 10, 2};
    fit->SetParameters(fit_p);
    fit->SetParNames("Norm1","Mean1", "SD1", "Norm2","Mean2", "SD2", "Norm3","Mean3", "SD3");
    fit->SetParLimits(0,1,1e10);
    fit->SetParLimits(3,1,1e10);
    fit->SetParLimits(6,1,1e10);
    fit->SetParLimits(1,-20,0);
    fit->SetParLimits(4,-10,10);
    fit->SetParLimits(7,0,20);
    fit->SetParLimits(2,0.1,10);
    fit->SetParLimits(5,0.1,10);
    fit->SetParLimits(8,0.1,10);    
    //fit->SetParLimits(9,0,1e10);    
    TFitResultPtr fitResult =h_vertex_proj->Fit(fit,"SIWL");
    
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
    
    //TF1* c = new TF1("const", "[0]",-20,20);
    //c->SetParameter(0,fitResult->Parameter(9));

    //c->SetLineColor(kCyan);
    //h_vertex_proj->GetListOfFunctions()->Add(c);

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
    double totalCounts =0; 
    double signalCounts=0;
    for (int i=0; i<nFixedz; i++){
        for (int j=0; j<nFixedz; j++){
            int binX1 = h_vertex_gem->GetXaxis()->FindBin(z_range[i]);
            int binX2 = h_vertex_gem->GetXaxis()->FindBin(z_range[i+1]);
            int binY1 = h_vertex_gem->GetYaxis()->FindBin(z_range[j]);
            int binY2 = h_vertex_gem->GetYaxis()->FindBin(z_range[j+1]);
            double content = h_vertex_gem->Integral(binX1, binX2, binY1, binY2);
            totalCounts+=content;
            if (i==j){
                signalCounts+=content;
            }
            std::cout << "Signal counts in region (" << z_range[i] << ", " << z_range[i+1] << ") x (" << z_range[j] << ", " << z_range[j+1] << "): " << content << std::endl;
        }
        TBox* box = new TBox(z_range[i],z_range[i], z_range[i+1], z_range[i+1]);
        box->SetLineColor(kRed);
        box->SetLineWidth(2);
        box->SetFillStyle(0);
        h_vertex_gem->GetListOfFunctions()->Add(box);
    }
    std::cout << "Total counts in vertex distribution: " << totalCounts << std::endl;
    std::cout << "Estimated signal counts in defined regions: " << signalCounts << std::endl;
    std::cout << "Estimated background counts: " << totalCounts - signalCounts << std::endl;
    std::cout << "Estimated signal-to-background ratio: " << signalCounts / (totalCounts - signalCounts) << std::endl;

    h_vertex_z->Fit(fit,"SIWL");
    


    //clean up
    outputFile->Write();
    outputFile->Close();
  
    std::cout << "Fitting completed. Results saved to " << outputFileName << std::endl;
}
