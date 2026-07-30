// -*- C++ -*-

#ifndef ANA_MANAGER_HH
#define ANA_MANAGER_HH

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>

#include <G4LorentzVector.hh>
#include <G4ThreeVector.hh>

#include <TParticle.h>
#include <TVector3.h>

class G4ParticleDefinition;

class TFile;
class TTree;

struct Track;

class VHitInfo;

//_____________________________________________________________________________
static const G4int MaxHits    = 500;
static const G4int MaxHitsTPC = 500;
static const G4int MaxPrimaryParticle = 10;

//const int MaxTrack = 1560*4;
//const int MaxTrack = 78*4;
const G4int MaxTrack = 54*20;

void initTrack(Track* aTrack);
void initTrack_ku(Track* aTrack);
int setInitialPara(Track* aTrack, double* initPara);
int setVirtualPlane(Track* aTrack);
void minuitInit(double printLevel);

static const int MAXtpctrNum=30;
static const int MAXtpctrhitNum=500;

//_____________________________________________________________________________
struct CounterData
{
  G4int ntrk;
  G4double resoX;
  G4int trackID;
  G4int particleID;
  G4int charge;
  G4double time;
  G4double beta;
  G4double edep;
  G4double dedx;
  G4double slength;
  G4double tlength;
  G4double mass;
  G4double pos0[3];
  G4double pos[3];
  G4double mom[4];
  G4int iLay;
  G4int iPad;
  G4int iRow;
  G4int parentID; //parent particle track id
  G4int parentPID;
  G4double res[3];
};

//_____________________________________________________________________________
struct TPCData
{
  G4int tpctr;
  G4int tpcpid;
  G4int tpcparentid;
  G4int tpcparentid_pid;
  G4double tpcpx;
  G4double tpcpy;
  G4double tpcpz;
  G4double tpcpp;

  G4double tpcpxfit;
  G4double tpcpyfit;
  G4double tpcpzfit;
  G4double tpcppfit;
  G4double tpcptfit;

  G4int tpcqq;
  G4double tpcpm;
  G4double tpcde;
  G4double tpclen;
  G4double tpcdedx;
  G4int tpclay;
  G4double tpcvtxpx;
  G4double tpcvtxpy;
  G4double tpcvtxpz;
  G4double tpcvtxx;
  G4double tpcvtxy;
  G4double tpcvtxz;
  //  G4double tpcene2
};

//_____________________________________________________________________________
struct Event
{
  Int_t evnum; // Event number
  Int_t generator;
  Int_t mode;        // mode number
  Int_t inc;        // INC id number

  Int_t HitNum_K;
  //  int tpctrNum_K;

  Int_t HitNum_p;
  //  int tpctrNum_p;

  std::map<TString, std::vector<TParticle>> hits;

  // TPCPad hit branches
  Int_t nhittpc;                    // Number of TPC pad hits stored in this event

  // Track and particle IDs
  std::vector<Int_t> ntrk;          // Track multiplicity/index value passed from the TPCPad hit
  std::vector<Int_t> trackidtpc;    // Geant4 track ID
  std::vector<Int_t> pidtpc;        // Particle PDG ID
  std::vector<Int_t> charge_tpc_hit; // Internal: Geant4 charge aligned with accepted pad hits
  std::vector<Int_t> parentidtpc;   // Parent Geant4 track ID
  std::vector<Int_t> parentpidtpc;  // Parent particle PDG ID

  // Hit positions
  std::vector<Double_t> xtpc;       // Smeared hit x position [mm]
  std::vector<Double_t> ytpc;       // Smeared hit y position [mm]
  std::vector<Double_t> ztpc;       // Smeared hit z position [mm]
  std::vector<Double_t> x0tpc;      // True Geant4 hit x position before smearing [mm]
  std::vector<Double_t> y0tpc;      // True Geant4 hit y position before smearing [mm]
  std::vector<Double_t> z0tpc;      // True Geant4 hit z position before smearing [mm]
  std::vector<Double_t> xtpc_pad;   // Pad-center x position from padtpc [mm]
  std::vector<Double_t> ytpc_pad;   // Pad-center y placeholder [mm]; currently copied from ytpc
  std::vector<Double_t> ztpc_pad;   // Pad-center z position from padtpc [mm]
  std::vector<Double_t> dxtpc_pad;  // True hit x minus pad-center x: x0tpc - xtpc_pad [mm]
  std::vector<Double_t> dytpc_pad;  // True hit y minus ytpc_pad: y0tpc - ytpc_pad; ytpc_pad currently equals ytpc [mm]
  std::vector<Double_t> dztpc_pad;  // True hit z minus pad-center z: z0tpc - ztpc_pad [mm]

  // Position resolution parameters
  std::vector<Double_t> resoX;      // Resolution scale used for smearing: transverse sigma in Simple mode, |ResVector| in Exp mode [mm]
  std::vector<Double_t> resxtpc;    // X component of the position resolution from GetResVector, not xtpc - x0tpc [mm]
  std::vector<Double_t> resytpc;    // Y component of the position resolution from GetResVector, not ytpc - y0tpc [mm]
  std::vector<Double_t> resztpc;    // Z component of the position resolution from GetResVector, not ztpc - z0tpc [mm]

  // Momentum and timing
  std::vector<Double_t> pxtpc;      // Particle momentum x at the hit [GeV/c]
  std::vector<Double_t> pytpc;      // Particle momentum y at the hit [GeV/c]
  std::vector<Double_t> pztpc;      // Particle momentum z at the hit [GeV/c]
  std::vector<Double_t> pptpc;      // Particle total momentum at the hit [GeV/c]
  std::vector<Double_t> timetpc;    // Global hit time [ns]
  std::vector<Double_t> tlengthtpc; // Geant4 track length at the hit [mm]
  std::vector<Double_t> betatpc;    // Particle beta at the hit

  // Energy loss and pad IDs
  std::vector<Double_t> edeptpc;    // Energy deposit in the pad (same as cluster de)
  std::vector<Double_t> dedxtpc;    // Energy deposit per pad path length
  std::vector<Double_t> slengthtpc; // Geant4 step length in the pad hit [mm]
  std::vector<Int_t> padtpc;       // Global pad ID calculated from layertpc and rowtpc
  std::vector<Int_t> layertpc;    // Pad layer ID
  std::vector<Int_t> rowtpc;      // Pad row ID

  Int_t ntrack;
  std::vector<Int_t> trackid;
  std::vector<Int_t> trackpid;
  std::vector<Int_t> charge_track;
  std::vector<Double_t> mom0_track;
  std::vector<Double_t> dedx_track;
  std::vector<Int_t> trackparentid;

  Int_t nvtx;
  std::vector<Int_t> vtx_type;       // 0: reaction, 1: decay
  std::vector<Int_t> vtx_motherpid;   // reaction: generator ID, decay: mother PDG
  std::vector<Double_t> vtx_x;
  std::vector<Double_t> vtx_y;
  std::vector<Double_t> vtx_z;
  std::vector<std::vector<Int_t>> vtx_trackid;
  std::vector<std::vector<Int_t>> vtx_trackpid;

  void ClearTPCPadHits()
  {
    nhittpc = 0;
    ntrk.clear();
    trackidtpc.clear();
    pidtpc.clear();
    charge_tpc_hit.clear();
    xtpc.clear();
    ytpc.clear();
    ztpc.clear();
    xtpc_pad.clear();
    ytpc_pad.clear();
    ztpc_pad.clear();
    dxtpc_pad.clear();
    dytpc_pad.clear();
    dztpc_pad.clear();
    x0tpc.clear();
    y0tpc.clear();
    z0tpc.clear();
    resoX.clear();
    resxtpc.clear();
    resytpc.clear();
    resztpc.clear();
    pxtpc.clear();
    pytpc.clear();
    pztpc.clear();
    pptpc.clear();
    timetpc.clear();
    tlengthtpc.clear();
    betatpc.clear();
    edeptpc.clear();
    dedxtpc.clear();
    slengthtpc.clear();
    padtpc.clear();
    layertpc.clear();
    rowtpc.clear();
    parentidtpc.clear();
    parentpidtpc.clear();
    ntrack = 0;
    trackid.clear();
    trackpid.clear();
    charge_track.clear();
    mom0_track.clear();
    dedx_track.clear();
    trackparentid.clear();
    nvtx = 0;
    vtx_type.clear();
    vtx_motherpid.clear();
    vtx_x.clear();
    vtx_y.clear();
    vtx_z.clear();
    vtx_trackid.clear();
    vtx_trackpid.clear();
  }

  void ResizeTPCPadHits(Int_t nhits)
  {
    nhittpc = nhits;
    ntrk.resize(nhits);
    trackidtpc.resize(nhits);
    pidtpc.resize(nhits);
    charge_tpc_hit.resize(nhits);
    xtpc.resize(nhits);
    ytpc.resize(nhits);
    ztpc.resize(nhits);
    xtpc_pad.resize(nhits);
    ytpc_pad.resize(nhits);
    ztpc_pad.resize(nhits);
    dxtpc_pad.resize(nhits);
    dytpc_pad.resize(nhits);
    dztpc_pad.resize(nhits);
    x0tpc.resize(nhits);
    y0tpc.resize(nhits);
    z0tpc.resize(nhits);
    resoX.resize(nhits);
    resxtpc.resize(nhits);
    resytpc.resize(nhits);
    resztpc.resize(nhits);
    pxtpc.resize(nhits);
    pytpc.resize(nhits);
    pztpc.resize(nhits);
    pptpc.resize(nhits);
    timetpc.resize(nhits);
    tlengthtpc.resize(nhits);
    betatpc.resize(nhits);
    edeptpc.resize(nhits);
    dedxtpc.resize(nhits);
    slengthtpc.resize(nhits);
    padtpc.resize(nhits);
    layertpc.resize(nhits);
    rowtpc.resize(nhits);
    parentidtpc.resize(nhits);
    parentpidtpc.resize(nhits);
  }
};

//_____________________________________________________________________________
class AnaManager
{
public:
  static G4String ClassName();
  static AnaManager& GetInstance();
  ~AnaManager();

private:
  AnaManager();
  AnaManager(const AnaManager&);
  AnaManager& operator=(const AnaManager&);

private:
  TFile* m_file;
  TTree* m_tree;
  TTree* m_tree_light;
  G4int m_on_off_helm;
  G4int m_pad_config;
  G4int m_experiment;
  G4bool m_do_accep_study;

  G4double m_effective_thickness;
  G4double m_mom_kaon_lab;
  G4double m_cos_theta;
  G4double m_cos_theta_lambda;
  
  CounterData counterData[MaxTrack];
  TPCData tpcData[MAXtpctrNum];

  int HitNum;
  int tpctrNum;
  int HitNum_K;
  int HitNum_p;

  G4double mean[MAXtpctrNum];//read fit parameters
  G4double trmean[MAXtpctrNum];//read fit parameters
  G4double cir_r[MAXtpctrNum];//read fit parameters
  G4double error[MAXtpctrNum];//read fit parameters
  G4double chi2[MAXtpctrNum];//read fit parameters
  G4double ndf[MAXtpctrNum];//read fit parameters
  G4double Pz[MAXtpctrNum];//read fit parameters

  ////////////////////getenv parameters
  G4double pad_length_in;
  G4double pad_length_out;
  G4double pad_gap;
  G4double pad_in_width;
  G4double pad_out_width;
  G4double pad_in_num;
  G4double pad_out_num;
  G4double truncated_mean_cut;

  G4double angle[40];
  G4double seg_angle[40];
  G4double seg_width[40];
  G4int numpads[40];
  G4double pad_in[40];
  G4double pad_out[40];
  G4double tpc_rad;

  // --------------------------------
  // combine beam and reaction generator  
  G4bool m_do_hit_tgt;
  G4bool m_do_generate_beam;
  G4bool m_do_combine;
  G4bool m_require_tpc_mp;
  G4bool m_threshold_con;
  G4int m_effective_evnum;
  G4int m_next_generator;
  G4int m_first_generator;
  G4int m_second_generator;
  G4ThreeVector m_next_pos;
  G4ThreeVector m_next_mom;
  G4ThreeVector m_vertex_pos;
  G4ThreeVector m_debug_pos;
  // --------------------------------

  // --------------------------------
  // for checking decay particle
  std::pair<G4String, G4String> m_previous_particle; // particle name, process name
  const std::unordered_map<G4int, G4String> m_focus_particle = {
    // generator id, particle name
    {7201, "kaon-"},
    {7202, "lambda"},
    {7203, "lambda"},
    {7204, "sigma-"},
    {7205, "lambda"},
    {7206, "sigma+"},
    {7207, "kaon-"},
    {7208, "kaon0S"},
    {7212, "lambda"},
    {7213, "lambda"},
    {7214, "sigma+"},
    {7215, "sigma-"}
  };
  G4int m_decay_particle_code;
  G4ThreeVector m_decay_position;
  // --------------------------------

  // --------------------------------
  // for acceptance study
  const G4double m_edep_threshold = 0.2; // MeV/cm
  const std::unordered_map<G4int, std::vector<G4int>> m_tpc_check_list = {
  //  gen   { check parentidtpc, PDG codes of check list }
    { 7202, {1, 2212, -211} }, // eta Lambda
    { 7203, {1, 2212, -211} }, // pi0 Lambda
    { 7204, {0, +211, -211} }, // pi+ Sigma-
    { 7205, {1, 2212, -211} }, // pi0 Sigma0
    { 7206, {0, -211, +211, 2212} }, // pi- Sigma+
    { 7207, {0, -321, 2212} }, // K p
    { 7208, {1, +211, -211} }, // k0 n
    { 7212, {0, 2212, -211, +211} }, // pi- pi+ Lambda
    { 7213, {0, 2212, -211, +211} }, // pi+ pi- Sigma0
    { 7214, {0, 2212, -211} }, // pi- pi0 Sigma+
    { 7215, {0, -211, +211} }, // pi+ pi0 Sigma-
    
    { 10402, {1, +321, -321, +321, -321}} //Phi Phi
  };
  const std::vector<G4int> m_forward_seg_narrow{16, 17, 18, 19, 20, 21, 22};
  const std::vector<G4int> m_forward_seg_wide{10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
  const std::vector<G4int> m_forward_seg_all{6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
  const G4double m_refractive_index_kvc = 1.46;
  const G4double m_refractive_index_bac = 1.115;

  G4int m_trig_flag_int;
  G4int m_focus_parent_id;  
  G4bool m_kaon_beam_flag;
  // --------------------------------

  
public:
  void BeginOfRunAction(G4int runnum);
  void EndOfRunAction();
  void BeginOfEventAction();
  int  EndOfEventAction();
  void MakeBranch(const G4String& sd_name);
  void MakeHistogram(const G4String& sd_name);
  void SetNhits(const G4String& sd_name, G4int nhits);
  void SetHitData(const VHitInfo* hit);
  void SetTPCData(G4int tpctr, G4int tpcpid, G4int tpcparentid,
                  G4int tpcparentid_pid, G4double tpcpx, G4double tpcpy,
                  G4double tpcpz,G4double tpcpp,  G4int tpcqq, G4double tpcpm,
                  G4double tpcde, G4double tpclen, G4int tpclay,
                  G4double vtxpxtpc2,G4double vtxpytpc2,G4double vtxpztpc2,
                  G4double vtxxtpc2,G4double vtxytpc2,G4double vtxztpc2,
                  G4double vtxenetpc2);
  void SetCounterDataSimple(G4int ntrk, G4double time, G4ThreeVector pos,
                      G4ThreeVector mom, G4int track, G4int particle,
                      G4int iLay, G4int iRow, G4double beta, G4double edep,
			    G4int parentid, G4int parentpid, G4int charge, G4double tlength, G4double slength);
  void SetCounterDataExp(G4int ntrk, G4double time, G4ThreeVector pos,
                      G4ThreeVector mom, G4int track, G4int particle,
                      G4int iLay, G4int iRow, G4double beta, G4double edep,
			 G4int parentid, G4int parentpid, G4int charge, G4double tlength, G4double slength);
  void SetFermiMomentum(const G4ThreeVector& p);
  void SetGeneratorID(G4int generator);
  void SetModeID(G4int mode);
  void SetIncID(G4int inc);
  void SetPrimaryParticle(G4int id, G4int pdg,
                          const G4LorentzVector& p,
                          const G4LorentzVector& v,
                          G4bool is_virtual_beam=false);
  void SetSecondaryVertex(G4int pdg, G4int motherPdg,
			  const G4LorentzVector& p,
			  const G4LorentzVector& v,
			  G4int daughterTrackID=0,
			  G4int motherTrackID=0);
  void SetBeamInfo(G4int pdg,
		   const G4LorentzVector& p,
		   const G4LorentzVector& v);
  void SetPrimaryVertex(G4int id, const G4ThreeVector& x);
  void SetPrimaryVertex(G4int id, G4double x, G4double y, G4double z);
  void SetEffectiveThickness(G4double effective_thickness);
  G4double GetEffectiveThickness();
  void SetMomKaonLab(G4double mom_kaon_lab);
  void SetCosTheta(G4double cos_theta);
  void SetCosThetaLambda(G4double cos_theta_lambda);
  void SetPreviousParticle(G4String particle_name, G4String process_name);
  
  // --------------------------------
  // combine beam and reaction generator
  void   SetDoHitTGT(G4bool do_hit_tgt);
  G4bool GetDoHitTGT();
  void   SetDoGenerateBeam(G4bool do_generate_beam);
  G4bool GetDoGenerateBeam();
  void   SetDoCombine(G4bool do_combine);
  G4bool GetDoCombine();
  void   SetRequireTpcMp(G4bool require_tpc_mp);
  G4bool GetRequireTpcMp();
  void   SetThresholdCondition(G4bool threshold_con);
  G4bool GetThresholdCondition();
  void  SetEffectiveEvnum(G4int effective_evnum);
  G4int GetEffectiveEvnum();
  void  SetNextGenerator(G4int next_generator);
  G4int GetNextGenerator();
  void  SetFirstGenerator(G4int first_generator);
  G4int GetFirstGenerator();
  void  SetSecondGenerator(G4int second_generator);
  G4int GetSecondGenerator();
  void          SetNextPos(G4double vx, G4double vy, G4double vz);
  G4ThreeVector GetNextPos();
  void          SetNextMom(G4double px, G4double py, G4double pz);
  G4ThreeVector GetNextMom();
  void          SetVertexPos(G4double vx, G4double vy, G4double vz);
  G4ThreeVector GetVertexPos();
  void          SetDebugPos(G4double vx, G4double vy, G4double vz);
  G4ThreeVector GetDebugPos();
  // --------------------------------

  // --------------------------------
  // for checking decay particle
  std::pair<G4String, G4String> GetPreviousParticle();
  G4String GetFocusParticle(G4int generator_id);
  void  SetDecayParticleCode(G4int decay_particle_code);
  G4int GetDecayParticleCode();
  void SetDecayPosition(G4ThreeVector decay_position);
  G4ThreeVector GetDecayPosition();
  G4bool IsInsideHtof(G4ThreeVector position);
  // --------------------------------

  // --------------------------------
  // for trigger
  void SetFocusParentID(G4int focus_parent_id);
  // --------------------------------

  void BuildVtxInfo();

  
  int CircleIntersect(double x1, double y1, double r1, double x2, double y2, double r2,
		      double ca1, double cb1, double ct01, int qq1,
		      double ca2, double cb2, double ct02, int qq2,
		      double inter1[3], double inter2[3])
  {
    // function inputs: x1, y1, r1, x2, y2, r2
    // function output: inter1, inter2 = coordinates of intersections

    double d,e,f,g,a,b,c;
    double x,y,discrim;

    if(x1 == x2 && y1 == y2){
      G4cout << x1 << " " << y1 << " " << r1 << G4endl;
      return 0;
    }
    //    G4cout << x1 << " " << y1 << " " << r1 << G4endl;
    //    G4cout << x2 << " " << y2 << " " << r2 << G4endl;

    d = -0.5*(r1*r1 - r2*r2 - x1*x1 - y1*y1 + x2*x2 + y2*y2);
    e =  0.5*(r1*r1 + r2*r2 - x1*x1 - y1*y1 - x2*x2 - y2*y2);
    if(fabs(y1-y2) < 1.0e-20) {
      x = d/(x1 - x2);
      a = 1.0;
      b = 0.0;
      c = x*x - x*(x1 + x2) - e;
      discrim = -4*a*c;
      if(discrim < 0) return 0;
      y = sqrt(discrim) / (2*a);
      inter1[0] = x;
      inter1[1] = y;
      inter2[0] = x;
      inter2[1] = -y;
      return 1;
    }
    f = (x1 - x2) / (y1 - y2);
    g = d / (y1 - y2);
    // cout << "d=" << d << " e=" << e << " f=" << f << " g=" << g << endl;

    a = 1. + f*f;
    b = f*(y1 + y2) - 2*f*g - (x1 + x2);
    c = g*g - g*(y1 + y2) - e;
    //    G4cout << "a=" << a << " b=" << b << " c=" << c << G4endl;

    discrim = b*b - 4*a*c;
    //    G4cout << "discrim = " << discrim << G4endl;
    if(discrim < 0) return 0;
    inter1[0] = (-b + sqrt(discrim)) / (2*a);
    inter1[1] = g - f*inter1[0];
    inter2[0] = (-b - sqrt(discrim)) / (2*a);
    inter2[1] = g - f*inter2[0];


    double theta11 = atan2(inter1[1]-y1, inter1[0]-x1);
    double theta21 = atan2(inter1[1]-y2, inter1[0]-x2);
    double theta12 = atan2(inter2[1]-y1, inter2[0]-x1);
    double theta22 = atan2(inter2[1]-y2, inter2[0]-x2);

    double tmp_y11 = -1.*(double)qq1*ca1*r1*(theta11-ct01)+cb1;
    double tmp_y21 = -1.*(double)qq2*ca2*r2*(theta21-ct02)+cb2;
    double tmp_y12 = -1.*(double)qq1*ca1*r1*(theta12-ct01)+cb1;
    double tmp_y22 = -1.*(double)qq2*ca2*r2*(theta22-ct02)+cb2;

    // std::cout<<"theta11="<<theta11<<", theta21="<<theta21
    // 	     <<", theta12="<<theta12<<", theta22="<<theta22<<std::endl;

    // std::cout<<"tmp_y11="<<tmp_y11<<", tmp_y21="<<tmp_y21<<std::endl;
    // std::cout<<"tmp_y12="<<tmp_y12<<", tmp_y22="<<tmp_y22<<std::endl;
    //getchar();

    inter1[2] = (tmp_y11+tmp_y21)/2.;
    inter2[2] = (tmp_y12+tmp_y22)/2.;


    return 1;
  }


  double linearFitter(const int np,
		      const double *x,
		      const double *y, double *er,
		      double *a, double *b){

    int i;

    double alpha=0.;
    double beta=0.;
    double gamma=0.;
    double AA=0;
    double BB=0;

    for(i=0;i<np;i++){
      alpha+=x[i]*x[i]/er[i]/er[i];
      beta+=x[i]/er[i]/er[i];
      gamma+=1./er[i]/er[i];
      AA+=y[i]*x[i]/er[i]/er[i];
      BB+=y[i]/er[i]/er[i];

      //  G4cout<<"x test: "<<x[i]<<G4endl;
      //  G4cout<<"y test: "<<y[i]<<G4endl;
    }

    //  G4cout<<"beta test: "<<beta<<G4endl;
    //  G4cout<<"alpha test: "<<alpha<<G4endl;
    //  G4cout<<"gamma test: "<<gamma<<G4endl;

    *a=(gamma*AA -  beta*BB)/(alpha*gamma-beta*beta);
    *b=(-beta *AA + alpha*BB)/(alpha*gamma-beta*beta);

    //  G4cout<<"a test: "<<(*a)<<G4endl;
    //  G4cout<<"b test: "<<(*b)<<G4endl;

    return 1.;
  }
};

//_____________________________________________________________________________
inline G4String
AnaManager::ClassName()
{
  static G4String s_name("AnaManager");
  return s_name;
}

//_____________________________________________________________________________
inline AnaManager&
AnaManager::GetInstance()
{
  static AnaManager s_instance;
  return s_instance;
}

#endif
