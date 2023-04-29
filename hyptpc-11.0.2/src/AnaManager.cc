// -*- C++ -*-

#include "AnaManager.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ParticleDefinition.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TParticle.h>
#include <TString.h>
#include <TTree.h>

#include "ConfMan.hh"
#include "DetectorConstruction.hh"
#include "FuncName.hh"
#include "ResHypTPC.hh"
#include "RungeKuttaTracker.hh"
#include "switch.h"
#include "track.hh"
#include "VHitInfo.hh"
#include "padHelper.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();
Event event;
std::map<TString, TH1*> hmap;
}

//_____________________________________________________________________________
AnaManager::AnaManager()
  : m_file(),
    m_tree(new TTree("g4hyptpc", "GEANT4 simulation for HypTPC"))
{
  m_tree->Branch("evnum", &event.evnum, "evnum/I");
  m_tree->Branch("generator", &event.generator, "generator/I");
  m_tree->Branch("mode",&event.mode,"mode/I");
  m_tree->Branch("inc",&event.inc,"inc/I");
  MakeBranch("PRM");
  for (const auto& sd_name : DetectorConstruction::GetSDList()) {
    MakeBranch(sd_name);
  }

  TString key;
  key = "Time";
  hmap[key] = new TH1D(key, key, 400, 0.0, 10.0);
  for(G4int i=0; i<G4ThreeVector::SIZE; ++i){
    key = Form("Pos%d", i);
    hmap[key] = new TH1D(key, key, 500, -25.0*CLHEP::cm, 25.0*CLHEP::cm);
    key = Form("Mom%d", i);
    if(i==2)
      hmap[key] = new TH1D(key, key, 400, 0.0, 2.0);
    else
      hmap[key] = new TH1D(key, key, 400, -2.0, 2.);
    key = Form("Fermi%d", i);
    hmap[key] = new TH1D(key, key, 500, -1.0*CLHEP::GeV, 1.0*CLHEP::GeV);
    hmap[key]->GetXaxis()->SetTitle("[MeV/c]");
  }

  return;
  m_tree->Branch("nhittpc",&event.nhittpc,"nhittpc/I");
  m_tree->Branch("ntrk",event.ntrk,"ntrk[nhittpc]/I");
  m_tree->Branch("ititpc",event.ititpc,"ititpc[nhittpc]/I");
  m_tree->Branch("idtpc",event.idtpc,"idtpc[nhittpc]/I");
  m_tree->Branch("xtpc",event.xtpc,"xtpc[nhittpc]/D");//after smeared by resolution
  m_tree->Branch("ytpc",event.ytpc,"ytpc[nhittpc]/D");//after smeared by resolution
  m_tree->Branch("ztpc",event.ztpc,"ztpc[nhittpc]/D");//after smeared by resolution
  m_tree->Branch("x0tpc",event.x0tpc,"x0tpc[nhittpc]/D");
  m_tree->Branch("y0tpc",event.y0tpc,"y0tpc[nhittpc]/D");
  m_tree->Branch("z0tpc",event.z0tpc,"z0tpc[nhittpc]/D");
  m_tree->Branch("resoX",event.resoX,"resoX[nhittpc]/D");
  m_tree->Branch("pxtpc",event.pxtpc,"pxtpc[nhittpc]/D");
  m_tree->Branch("pytpc",event.pytpc,"pytpc[nhittpc]/D");
  m_tree->Branch("pztpc",event.pztpc,"pztpc[nhittpc]/D");
  m_tree->Branch("pptpc",event.pptpc,"pptpc[nhittpc]/D");   // total mometum
  m_tree->Branch("masstpc",event.masstpc,"masstpc[nhittpc]/D");   // mass TPC
  m_tree->Branch("timetpc",event.timetpc,"timetpc[nhittpc]/D");
  m_tree->Branch("betatpc",event.betatpc,"betatpc[nhittpc]/D");
  m_tree->Branch("edeptpc",event.edeptpc,"edeptpc[nhittpc]/D");
  m_tree->Branch("dedxtpc",event.dedxtpc,"dedxtpc[nhittpc]/D");
  m_tree->Branch("slengthtpc",event.slengthtpc,"slengthtpc[nhittpc]/D");
  m_tree->Branch("tlengthtpc",event.tlengthtpc,"tlengthtpc[nhittpc]/D");
  m_tree->Branch("iPadtpc",event.iPadtpc,"iPadtpc[nhittpc]/I");
  m_tree->Branch("laytpc",event.laytpc,"laytpc[nhittpc]/I");
  m_tree->Branch("rowtpc",event.rowtpc,"rowtpc[nhittpc]/I");
  m_tree->Branch("parentID",event.parentID,"parentID[nhittpc]/I");
  m_tree->Branch("xtpc_pad",event.xtpc_pad,"xtpc_pad[nhittpc]/D");//pad center position
  m_tree->Branch("ytpc_pad",event.ytpc_pad,"ytpc_pad[nhittpc]/D");//pad center position (dummy = ytpc)
  m_tree->Branch("ztpc_pad",event.ztpc_pad,"ztpc_pad[nhittpc]/D");//pad center position
  m_tree->Branch("dxtpc_pad",event.dxtpc_pad,"dxtpc_pad[nhittpc]/D");//x0tpc - xtpc_pad
  m_tree->Branch("dytpc_pad",event.dytpc_pad,"dytpc_pad[nhittpc]/D");//y0tpc - ytpc_pad (dummy = 0)
  m_tree->Branch("dztpc_pad",event.dztpc_pad,"dztpc_pad[nhittpc]/D");//z0tpc - ztpc_pad

  //// Study on multiplicity
  // m_tree->Branch("nthlay",event.nthlay,"nthlay[nhittpc]/I");
  // m_tree->Branch("nthpad",event.nthpad,"nthpad[nhittpc]/I");
  // m_tree->Branch("laypad",event.laypad,"laytpadpc[nhittpc][40][250]/I");


  //shhwang ntrtpc --> number of trak in tpc
  // m_tree->Branch("ntrtpc",&event.ntrtpc,"ntrtpc/I");
  // m_tree->Branch("trpmtpc",event.trpmtpc,"trpmtpc[ntrtpc]/D");
  // m_tree->Branch("trqqtpc",event.trqqtpc,"trqqtpc[ntrtpc]/I");
  // m_tree->Branch("trpidtpc",event.trpidtpc,"trpidtpc[ntrtpc]/I");
  // m_tree->Branch("trparentidtpc",event.trparentidtpc,"trparentidtpc[ntrtpc]/I");
  // //m_tree->Branch("trparentid_pid_tpc",event.trparentid_pid_tpc,"trparentid_pid_tpc[ntrtpc]/I");

  // m_tree->Branch("trpxtpc",event.trpxtpc,"trpxtpc[ntrtpc]/D");
  // m_tree->Branch("trpytpc",event.trpytpc,"trpytpc[ntrtpc]/D");
  // m_tree->Branch("trpztpc",event.trpztpc,"trpztpc[ntrtpc]/D");
  // m_tree->Branch("trpptpc",event.trpptpc,"trpptpc[ntrtpc]/D");
  // m_tree->Branch("trpttpc",event.trpttpc,"trpttpc[ntrtpc]/D");

  // m_tree->Branch("trpxtpcfit",event.trpxtpcfit,"trpxtpcfit[ntrtpc]/D");
  // m_tree->Branch("trpytpcfit",event.trpytpcfit,"trpytpcfit[ntrtpc]/D");
  // m_tree->Branch("trpztpcfit",event.trpztpcfit,"trpztpcfit[ntrtpc]/D");
  // m_tree->Branch("trpptpcfit",event.trpptpcfit,"trpptpcfit[ntrtpc]/D");
  // m_tree->Branch("trpttpcfit",event.trpttpcfit,"trpttpcfit[ntrtpc]/D");

  // m_tree->Branch("vtpxtpc",event.vtpxtpc,"vtpxtpc[ntrtpc]/D");
  // m_tree->Branch("vtpytpc",event.vtpytpc,"vtpytpc[ntrtpc]/D");
  // m_tree->Branch("vtpztpc",event.vtpztpc,"vtpztpc[ntrtpc]/D");
  // m_tree->Branch("vtpptpc",event.vtpptpc,"vtpptpc[ntrtpc]/D");

  // m_tree->Branch("vtxtpc",event.vtxtpc,"vtxtpc[ntrtpc]/D");
  // m_tree->Branch("vtytpc",event.vtytpc,"vtytpc[ntrtpc]/D");
  // m_tree->Branch("vtztpc",event.vtztpc,"vtztpc[ntrtpc]/D");

  // m_tree->Branch("vtxtpcfit",event.vtxtpcfit,"vtxtpcfit[ntrtpc]/D");
  // m_tree->Branch("vtytpcfit",event.vtytpcfit,"vtytpcfit[ntrtpc]/D");
  // m_tree->Branch("vtztpcfit",event.vtztpcfit,"vtztpcfit[ntrtpc]/D");

  // m_tree->Branch("trdetpc",event.trdetpc,"trdetpc[ntrtpc]/D");
  // m_tree->Branch("trlentpc",event.trlentpc,"trlentpc[ntrtpc]/D");
  // m_tree->Branch("trdedxtpc",event.trdedxtpc,"trdedxtpc[ntrtpc]/D");
  // m_tree->Branch("trdedxtrtpc",event.trdedxtrtpc,"trdedxtrtpc[ntrtpc]/D");
  // m_tree->Branch("trlaytpc",event.trlaytpc,"trlaytpc[ntrtpc]/I");
  // m_tree->Branch("cir_r",event.cir_r,"cir_r[ntrtpc]/D");
  // m_tree->Branch("cir_x",event.cir_x,"cir_x[ntrtpc]/D");
  // m_tree->Branch("cir_z",event.cir_z,"cir_z[ntrtpc]/D");
  // m_tree->Branch("cir_fit",event.cir_fit,"cir_fit[ntrtpc]/D");
  // m_tree->Branch("vtx_flag",event.vtx_flag,"vtx_flag[ntrtpc]/I");
  // m_tree->Branch("a_fory",event.a_fory,"a_fory[ntrtpc]/D");
  // m_tree->Branch("b_fory",event.b_fory,"b_fory[ntrtpc]/D");
}

//_____________________________________________________________________________
AnaManager::~AnaManager()
{
  if (gFile && gFile->IsOpen())
    gFile->Close();
}

//_____________________________________________________________________________
void
AnaManager::MakeBranch(const G4String& sd_name)
{
  static const Int_t bufsize = 32000;
  m_tree->Branch(sd_name.data(),
               "std::vector<TParticle>",
               &event.hits[sd_name], bufsize, -1);
}

//_____________________________________________________________________________
void
AnaManager::BeginOfRunAction(G4int /* runnum */)
{
  m_file = new TFile(gConf.Get<G4String>("ROOT"), "RECREATE");
  static auto obj = new TNamed("conf", gConf.ConfBuf());
  obj->Write();
  m_tree->Reset();

  for (auto& h: hmap) {
    h.second->Reset();
  }

  event.evnum = 0;

  G4double target_pos_z=-143.;
  truncated_mean_cut = gConf.Get<G4double>("TruncatedMeanCut");
  m_experiment = gConf.Get<G4int>("Experiment");
  //out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  pad_length_in = gConf.Get<G4double>("PadLengthIn");
  pad_length_out = gConf.Get<G4double>("PadLengthOut");
  pad_gap = gConf.Get<G4double>("PadGap");

  ////pad configure
  m_pad_config = gConf.Get<G4int>("PadConfigure");
  pad_in_num = gConf.Get<G4int>("PadNumIn");
  pad_out_num = gConf.Get<G4int>("PadNumOut");
  pad_in_width = gConf.Get<G4double>("PadWidthOut");
  pad_out_width = gConf.Get<G4double>("PadWidthOut");

  m_on_off_helm = gConf.Get<G4int>("ShsFieldMap");

  for(G4int i=0.;i<40;i++){
    angle[i]=0;
    seg_angle[i]=0;
    seg_width[i]=0;
    numpads[i]=0;

    pad_in[i]=0;
    pad_out[i]=0;
  }
  tpc_rad=250;
  G4double cen_diff=fabs(target_pos_z);

  if(m_pad_config ==1){
    for(G4int i=0;i<pad_in_num+pad_out_num;i++){
      if(i<pad_in_num){
	pad_in[i]=10.+(pad_length_in+pad_gap)*i;
	pad_out[i]=10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]=360.;
      }else {
	pad_in[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num);
	pad_out[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num) + pad_length_out;
	angle[i]=180.-acos((pow(pad_out[i],2)+pow(cen_diff,2)-pow(tpc_rad,2))/(2*pad_out[i]*cen_diff))*180./acos(-1.);
      }
      //      G4cout<<angle[i]<<G4endl;
      //      G4cout<<pad_in[i]<<G4endl;
    }


  }else if(m_pad_config ==2){
    for(G4int i=0;i<pad_in_num+pad_out_num;i++){
      if(i<pad_in_num){
	pad_in[i]=10.+(pad_length_in+pad_gap)*i;
	pad_out[i]=10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]=360.;
	if(i==0){
	  numpads[i]=48.;
	}else if(i<pad_in_num){
	  numpads[i]=24.*2.*(i+1.)/2.;
	}
      }else {
	pad_in[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num);
	pad_out[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num) + pad_length_out;
      }
    }
    angle[10]=180.-155.35;
    angle[11]=180.-144.8;
    angle[12]=180.-138.;
    angle[13]=180.-116.73;
    angle[14]=180.-106.;
    angle[15]=180.-98.77;
    angle[16]=180.-94.29;
    angle[17]=180.-89.8;
    angle[18]=180.-87.18;
    angle[19]=180.-84.16;
    angle[20]=180.-81.48;
    angle[21]=180.-73.39;
    angle[22]=180.-65.51011;
    angle[23]=180.-60.19;
    angle[24]=180.-56.35239;
    angle[25]=180.-52.85;
    angle[26]=180.-50.14;
    angle[27]=180.-47.17;
    angle[28]=180.-41.24;
    angle[29]=180.-29.;
    angle[30]=180.-23.23;
    angle[31]=180.-18.69;

    numpads[10]=208.;
    numpads[11]=218.;
    numpads[12]=230.;
    numpads[13]=214.;
    numpads[14]=212.;
    numpads[15]=214.;
    numpads[16]=220.;
    numpads[17]=224.;
    numpads[18]=232.;
    numpads[19]=238.;
    numpads[20]=244.;
    numpads[21]=232.;
    numpads[22]=218.;
    numpads[23]=210.;
    numpads[24]=206.;
    numpads[25]=202.;
    numpads[26]=200.;
    numpads[27]=196.;
    numpads[28]=178.;
    numpads[29]=130.;
    numpads[30]=108.;
    numpads[31]=90.;
    G4int all_channels=0;
    G4int all_channels2=0;
    G4int num_pad_check=0;


    for(G4int i=0;i<pad_in_num+pad_out_num;i++){
      if(i<pad_in_num){
	seg_angle[i]=360./double(numpads[i]);
	seg_width[i]=pad_in[i]*(angle[i])*CLHEP::pi/180./numpads[i];

	num_pad_check=angle[i]/seg_angle[i];
      }else if(i>=pad_in_num){
	seg_angle[i]=(180.-angle[i])*2/double(numpads[i]);
	seg_width[i]=pad_in[i]*(180-angle[i])*2.*acos(-1.)/180./numpads[i];
	num_pad_check=(180.-angle[i])*2/seg_angle[i];
      }

      G4cout<<i<<" degree :"<<seg_angle[i]<<G4endl;
      G4cout<<i<<" width :"<<seg_angle[i]*acos(-1.)/180.*pad_in[i]<<G4endl;

      all_channels=all_channels+numpads[i];
      all_channels2=all_channels2+num_pad_check;
    }
    G4cout<<"------------------------"<<G4endl;
    G4cout<<"Total pads:"<<all_channels<<G4endl;
    G4cout<<"Total pads(check):"<<all_channels<<G4endl;
    G4cout<<"------------------------"<<G4endl;
  }
}

//_____________________________________________________________________________
void
AnaManager::EndOfRunAction()
{
  m_file->cd();
  m_tree->Write();
  for (auto& h: hmap) {
    h.second->Write();
  }
  m_file->Close();
}

//_____________________________________________________________________________
void
AnaManager::BeginOfEventAction()
{
  HitNum=0;
  tpctrNum=0;


  //for K+
  HitNum_K=0;
  //  tpctrNum_K=0;

  //for proton
  HitNum_p=0;
  //  tpctrNum_K=0;

  event.hits["PRM"].clear();
  for (const auto& sd_name: DetectorConstruction::GetSDList()) {
    event.hits[sd_name].clear();
  }

  event.nhittpc = 0;
  event.ntrtpc = 0;

  event.HitNum_K=-1;

  event.HitNum_p=-1;


  /* ntrtpc initialization */

  for(G4int i=0; i<MaxHitsTPC;++i){
    event.trpidtpc[i]  = -1;
    event.trparentidtpc[i]  = -1;
    event.trparentid_pid_tpc[i]  = -1;


    event.trpptpc[i]  = -9999.9999;
    event.trpttpc[i]  = -9999.9999;
    event.trpxtpc[i]  = -9999.9999;
    event.trpytpc[i]  = -9999.9999;
    event.trpztpc[i]  = -9999.9999;

    event.vtpxtpc[i]  = -9999.9999;
    event.vtpytpc[i]  = -9999.9999;
    event.vtpztpc[i]  = -9999.9999;

    event.vtxtpc[i]  = -9999.9999;
    event.vtytpc[i]  = -9999.9999;
    event.vtztpc[i]  = -9999.9999;

    event.trpttpcfit[i]  = -9999.9999;

    event.trpptpcfit[i]  = -9999.9999;
    event.trpxtpcfit[i]  = -9999.9999;
    event.trpytpcfit[i]  = -9999.9999;
    event.trpztpcfit[i]  = -9999.9999;

    event.trpmtpc[i]  = -9999.9999;
    event.trqqtpc[i]  = -9999;

    event.trdetpc[i]  = -9999.9999;
    event.trlentpc[i]  = -9999.9999;
    event.trdedxtpc[i]  = -9999.9999;
    event.trdedxtrtpc[i]  = -9999.9999;
    event.trlaytpc[i]  = -9999;

    event.cir_r[i]  = -9999.9999;
    event.cir_x[i]  = -9999.9999;
    event.cir_z[i]  = -9999.9999;
    event.cir_fit[i]  = -9999.9999;

    event.vtx_flag[i]  = -1;
    event.a_fory[i]  = -9999.9999;
    event.b_fory[i]  = -9999.9999;
  }


  for(int i=0;i<MaxTrack;i++){

    /// initialization pad multiplicity
    event.nthlay[i]=-9999.;
    event.nthpad[i]=-9999.;
    for(int j = 0; j< MaxNthLay;j++){
      for(int k = 0; k< MaxNthPad;k++){
	event.laypad[i][j][k]  = 0.;
      }
    }
    //////////////


    event.xtpc[i] = -9999.9;
    event.ytpc[i] = -9999.9;
    event.ztpc[i] = -9999.9;

    event.xtpc_pad[i] = -9999.9;
    event.ytpc_pad[i] = -9999.9;
    event.ztpc_pad[i] = -9999.9;
    event.dxtpc_pad[i] = -9999.9;
    event.dytpc_pad[i] = -9999.9;
    event.dztpc_pad[i] = -9999.9;



    event.x0tpc[i] = -9999.9;
    event.y0tpc[i] = -9999.9;
    event.z0tpc[i] = -9999.9;
    event.resoX[i] = -9999.9;

    event.pxtpc[i] = -9999.9;
    event.pytpc[i] = -9999.9;
    event.pztpc[i] = -9999.9;
    event.pptpc[i] = -9999.9;

    event.masstpc[i] = -9999.9;

    event.timetpc[i] = -9999.9;
    event.betatpc[i] = -9999.9;

    event.edeptpc[i] = -9999.9;

    event.ititpc[i] = -1;
    event.idtpc[i] = -1;
    event.iPadtpc[i] = -1;
    event.laytpc[i] = -1;
    event.rowtpc[i] = -1;
    event.parentID[i] = -1;
  }
}

//_____________________________________________________________________________
int
AnaManager::EndOfEventAction()
{
  event.evnum++;

  if(HitNum > 0){
    G4int c[MAX_TRACK] = {};

    for(G4int i=0;i<MAX_TRACK;i++){
      mean[i]=0.;
      trmean[i]=0.;
    }

    G4double vtxxfit[MAX_TRACK];//read fit parameters
    G4double vtxyfit[MAX_TRACK];//read fit parameters
    G4double vtxzfit[MAX_TRACK];//read fit parameters

    G4double vtxpxfit[MAX_TRACK];//read fit parameters
    // G4double vtxpyfit[MAX_TRACK];//read fit parameters
    G4double vtxpzfit[MAX_TRACK];//read fit parameters

    for(G4int i=0;i<MAX_TRACK;i++){
      vtxxfit[i]=-9999.9999;
      vtxyfit[i]=-9999.9999;
      vtxzfit[i]=-9999.9999;
      vtxpxfit[i]=-9999.9999;
      // vtxpyfit[i]=-9999.9999;
      vtxpzfit[i]=-9999.9999;
      Pz[i]=-9999.9999;
    }

    G4double x[MAX_TRACK][MAXtpctrhitNum]={{-9999.9999},{-9999.9999}};;
    G4double z[MAX_TRACK][MAXtpctrhitNum]={{-9999.9999},{-9999.9999}};
    G4double y[MAX_TRACK][MAXtpctrhitNum]={{-9999.9999},{-9999.9999}};
    G4double ede[MAX_TRACK][MAXtpctrhitNum]={{0.},{0.}};

    ////// shhwang position read
    ///shhwang code


    if(tpctrNum>9){
      G4cout<<"Error--> over the number of tracks in the TPC:"<<tpctrNum<<G4endl;
    }

    G4int sh_paID[MAX_TRACK] = {};
    for(G4int i=0; i<HitNum; i++){
      G4int ii=counterData[i].ntrk;
      x[ii][c[ii]]=counterData[i].pos[0];
      z[ii][c[ii]]=counterData[i].pos[2];
      y[ii][c[ii]]=counterData[i].pos[1];
      sh_paID[ii]=counterData[i].parentID;
      ede[ii][c[ii]]=counterData[i].dedx;
      c[ii]=c[ii]+1;
    }

    G4double test[MAX_TRACK]={-1};
    G4double cx[MAX_TRACK]={-9999.9999};
    G4double cz[MAX_TRACK]={-9999.9999};
    G4double cir_x[MAX_TRACK]={-9999.9999};
    G4double cir_z[MAX_TRACK]={-9999.9999};
    G4double rad[MAX_TRACK]={-9999.9999};
    G4double a_fory[MAX_TRACK]={-9999.9999};
    G4double b_fory[MAX_TRACK]={-9999.9999};
    G4double theta0_fory[MAX_TRACK]={-9999.9999};
    G4int vtx_flag[MAX_TRACK]={-1};

    for(G4int i=0;i<MAX_TRACK;i++){
      cir_r[i]=-9999.9999;
      cir_x[i]=-9999.9999;
      cir_z[i]=-9999.9999;
      mean[i]=-9999.9999;
    }



    for(G4int kk=0; kk<tpctrNum; kk++){
      if(c[kk]>3.){
	//	G4cout<<"start circle fit"<<G4endl;
	// test[kk]=circleFit(x[kk],z[kk],y[kk],c[kk],&cx[kk],&cz[kk],&rad[kk],&Pz[kk],
	// 		   &a_fory[kk], &b_fory[kk], &theta0_fory[kk]);
	if(test[kk]!=-1.){
	  cir_r[kk]=rad[kk];
	  cir_x[kk]=cx[kk];
	  cir_z[kk]=cz[kk];
	}

      }
    }
    G4double mom_theta[MAX_TRACK]={0.};

    // calcute vtx with production points
    for(G4int i=0;i<MAX_TRACK;i++){
      G4double rho1 = rad[i];
      G4double cx1 = cx[i];
      G4double cz1 = cz[i];
      G4double cx2 = event.hits["PRM"][0].Vx();
      G4double cz2 = event.hits["PRM"][0].Vz();
      G4double theta12=atan2(cz2-cz1, cx2-cx1);
      G4double ca1=a_fory[i];
      G4double cb1=b_fory[i];
      G4double ct01=theta0_fory[i];


      G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

      vtxxfit[i]=cos(theta12)*cent_dist+cx1;
      vtxzfit[i]=sin(theta12)*cent_dist+cz1;
      vtxyfit[i]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

      mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;

      // vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
      // vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);

      vtx_flag[i]=1;
    }

    ////think about parent ID
    ////--> find the track with same parent ID
    //// sh_

    for(G4int i=0;i<MAX_TRACK;i++){
      for(G4int j=i;j<MAX_TRACK;j++){
	if(i!=j && (test[i]>0 && test[j]>0)){
	  if(sh_paID[i]==sh_paID[j] && sh_paID[i]>0. && sh_paID[j]>0.){
	    //	    G4cout<<"vtx1"<<env_helm_field<<G4endl;
	    G4double rho1=rad[i];
	    G4double rho2=rad[j];

	    G4double cx1=cx[i];
	    G4double cz1=cz[i];
	    G4double ca1=a_fory[i];
	    G4double cb1=b_fory[i];
	    G4double ct01=theta0_fory[i];

	    G4double cx2=cx[j];
	    G4double cz2=cz[j];
	    G4double ca2=a_fory[j];
	    G4double cb2=b_fory[j];
	    G4double ct02=theta0_fory[j];


	    G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

	    double point1[3]={0};
	    double point2[3]={0};
	    G4int k;

	    if((cent_dist-(rho1+rho2))>0.){


	      G4double theta12=atan2(cz2-cz1,cx2-cx1);
	      G4double centr=rho1+(cent_dist-(rho1+rho2))/2;

	      point1[0]=cos(theta12)*centr+cx1;
	      point1[1]=sin(theta12)*centr+cz1;
	      point1[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

	      G4double theta21=atan2(cz1-cz2,cx1-cx2);
	      G4double centr1=rho2+(cent_dist-(rho1+rho2))/2;
	      point2[0]=cos(theta21)*centr1+cx2;
	      point2[1]=sin(theta21)*centr1+cz2;
	      point2[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta21-ct02)+cb2;

	      vtxxfit[i]=point1[0];
	      vtxzfit[i]=point1[1];
	      vtxyfit[i]=(point1[2]+point2[2])/2.;
	      vtxxfit[j]=point2[0];
	      vtxzfit[j]=point2[1];
	      vtxyfit[j]=(point1[2]+point2[2])/2.;
	      vtx_flag[i]=2;
	      vtx_flag[j]=2;
	      //
	    }else  if((cent_dist+fmin(rho1,rho2))<fmax(rho1,rho2)){
	      if(rho1>=rho2){ //rho1>rho2
		G4double theta12=atan2(cz2-cz1,cx2-cx1);
		G4double centr=rho1-(rho1-cent_dist-rho2)/2; //rho1>rho2
		point1[0]=cos(theta12)*centr+cx1;
		point1[1]=sin(theta12)*centr+cz1;
		point1[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

		G4double theta21=atan2(cz2-cz1,cx2-cx1);
		G4double centr1=rho2+(rho1-cent_dist-rho2)/2.; //rho1>rho2
		point2[0]=cos(theta21)*centr1+cx2;
		point2[1]=sin(theta21)*centr1+cz2;
		point2[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta21-ct02)+cb2;
		//		G4cout<<"test1"<<G4endl;

	      }else if(rho2>rho1){ //rho1<rho2
		G4double theta12=atan2(cz1-cz2,cx1-cx2);
		G4double centr=rho2-(rho2-cent_dist-rho1)/2; //rho1<rho2
		point1[0]=cos(theta12)*centr+cx2;
		point1[1]=sin(theta12)*centr+cz2;
		point1[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta12-ct02)+cb2;

		G4double theta21=atan2(cz1-cz2,cx1-cx2);
		G4double centr1=rho1+(rho2-cent_dist-rho1)/2; //rho1<rho2
		point2[0]=cos(theta21)*centr1+cx1;
		point2[1]=sin(theta21)*centr1+cz1;
		point2[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta21-ct01)+cb1;
	      }

	      vtxxfit[i]=point1[0];
	      vtxzfit[i]=point1[1];
	      vtxyfit[i]=(point1[2]+point2[2])/2.;
	      // vtxxfit[j]=point1[0];
	      // vtxzfit[j]=point1[1];
	      // vtxyfit[j]=point1[2];
	      vtxxfit[j]=point2[0];
	      vtxzfit[j]=point2[1];
	      vtxyfit[j]=(point1[2]+point2[2])/2.;

	      vtx_flag[i]=3;
	      vtx_flag[j]=3;
	    } else {

	      //k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,point1,point2);
	      k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,ca1,cb1,ct01,tpcData[i].tpcqq,ca2,cb2,ct02,tpcData[j].tpcqq,point1,point2);
	      if(k == 0) {
		G4cout << "no solution" << G4endl;
	      }else if(k>0){


		G4double dist1=sqrt(pow(point1[0]-tpcData[i].tpcvtxx,2)+pow(point1[1]-tpcData[i].tpcvtxz,2));
		G4double dist2=sqrt(pow(point2[0]-tpcData[i].tpcvtxx,2)+pow(point2[1]-tpcData[i].tpcvtxz,2));

		if(dist1<=dist2){//point1 is correct
		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];
		  vtxyfit[i]=point1[2];

		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  vtxyfit[j]=point1[2];
		}else if(dist1>dist2){//point1 is correct
		  vtxxfit[i]=point2[0];
		  vtxzfit[i]=point2[1];
		  vtxyfit[i]=point2[2];

		  vtxxfit[j]=point2[0];
		  vtxzfit[j]=point2[1];
		  vtxyfit[j]=point2[2];
		}
		vtx_flag[i]=4;
		vtx_flag[j]=4;
	      }

	      mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;
	      mom_theta[j]=atan2(vtxzfit[j]-cz[j],vtxxfit[j]-cx[j])-acos(-1.)/2;


	      // std::cout<<"x01="<<x[i][0]<<", x2="<<x[j][0]<<std::endl;
	      // std::cout<<"y01="<<y[i][0]<<", y2="<<y[j][0]<<std::endl;
	      // std::cout<<"z01="<<z[i][0]<<", z2="<<z[j][0]<<std::endl;

	      // std::cout<<"vtx fit="<<vtxxfit[i]<<", true vtx="<<tpcData[i].tpcvtxx<<std::endl;
	      // std::cout<<"vty fit="<<vtxyfit[i]<<", true vty="<<tpcData[i].tpcvtxy<<std::endl;
	      // std::cout<<"vtz fit="<<vtxzfit[i]<<", true vtz="<<tpcData[i].tpcvtxz<<std::endl;
	      //getchar();


	      // vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      // vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      // vtxpxfit[j]=cos(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      // vtxpzfit[j]=sin(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      //	    G4cout<<"bfield:"<<env_helm_field<<G4endl;


	    }

	    ///from vertex particle, but it need more than 2
	  }else if(sh_paID[i]==sh_paID[j] && sh_paID[i]==0. && sh_paID[j]==0.){
	    //	    G4cout<<"vtx2"<<env_helm_field<<G4endl;
	    G4double rho1=rad[i];
	    G4double rho2=rad[j];

	    G4double cx1=cx[i];
	    G4double cz1=cz[i];
	    G4double ca1=a_fory[i];
	    G4double cb1=b_fory[i];
	    G4double ct01=theta0_fory[i];

	    G4double cx2=cx[j];
	    G4double cz2=cz[j];
	    G4double ca2=a_fory[j];
	    G4double cb2=b_fory[j];
	    G4double ct02=theta0_fory[j];

	    G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

	    double point1[3]={0};
	    double point2[3]={0};
	    G4int k;

	    if((cent_dist-(rho1+rho2))>0.){


	      G4double theta12=atan2(cz2-cz1,cx2-cx1);
	      G4double centr=rho1+(cent_dist-(rho1+rho2))/2;

	      point1[0]=cos(theta12)*centr+cx1;
	      point1[1]=sin(theta12)*centr+cz1;
	      point1[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

	      G4double theta21=atan2(cz1-cz2,cx1-cx2);
	      G4double centr1=rho2+(cent_dist-(rho1+rho2))/2;
	      point2[0]=cos(theta21)*centr1+cx2;
	      point2[1]=sin(theta21)*centr1+cz2;
	      point2[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta21-ct02)+cb2;

	      vtxxfit[i]=point1[0];
	      vtxzfit[i]=point1[1];
	      vtxyfit[i]=(point1[2]+point2[2])/2.;
	      vtxxfit[j]=point1[0];
	      vtxzfit[j]=point1[1];
	      vtxyfit[j]=(point1[2]+point2[2])/2.;

	      vtx_flag[i]=5;
	      vtx_flag[j]=5;
	    }else  if((cent_dist+fmin(rho1,rho2))<fmax(rho1,rho2)){

	      if(rho1>=rho2){ //rho1>rho2
		G4double theta12=atan2(cz2-cz1,cx2-cx1);
		G4double centr=rho1-(rho1-cent_dist-rho2)/2; //rho1>rho2
		point1[0]=cos(theta12)*centr+cx1;
		point1[1]=sin(theta12)*centr+cz1;
		point1[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

		G4double theta21=atan2(cz2-cz1,cx2-cx1);
		G4double centr1=rho2+(rho1-cent_dist-rho2)/2.; //rho1>rho2
		point2[0]=cos(theta21)*centr1+cx2;
		point2[1]=sin(theta21)*centr1+cz2;
		point2[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta21-ct02)+cb2;
		//		G4cout<<"test1"<<G4endl;

	      }else if(rho2>rho1){ //rho1<rho2
		G4double theta12=atan2(cz1-cz2,cx1-cx2);
		G4double centr=rho2-(rho2-cent_dist-rho1)/2; //rho1<rho2
		point1[0]=cos(theta12)*centr+cx2;
		point1[1]=sin(theta12)*centr+cz2;
		point1[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta12-ct02)+cb2;

		G4double theta21=atan2(cz1-cz2,cx1-cx2);
		G4double centr1=rho1+(rho2-cent_dist-rho1)/2; //rho1<rho2
		point2[0]=cos(theta21)*centr1+cx1;
		point2[1]=sin(theta21)*centr1+cz1;
		point2[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta21-ct01)+cb1;
	      }

	      vtxxfit[i]=point1[0];
	      vtxzfit[i]=point1[1];
	      vtxyfit[i]=(point1[2]+point2[2])/2.;
	      vtxxfit[j]=point2[0];
	      vtxzfit[j]=point2[1];
	      vtxyfit[j]=(point1[2]+point2[2])/2.;
	      vtx_flag[i]=6;
	      vtx_flag[j]=6;
	    } else {

	      //k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,point1,point2);
	      k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,ca1,cb1,ct01,tpcData[i].tpcqq,ca2,cb2,ct02,tpcData[j].tpcqq,point1,point2);
	      if(k == 0) {
		G4cout << "no solution" << G4endl;
	      }else if(k>0){

		G4double dist1=sqrt(pow(point1[0]-tpcData[i].tpcvtxx,2)+pow(point1[1]-tpcData[i].tpcvtxz,2));
		G4double dist2=sqrt(pow(point2[0]-tpcData[i].tpcvtxx,2)+pow(point2[1]-tpcData[i].tpcvtxz,2));

		if(dist1<=dist2){//point1 is correct
		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];
		  vtxyfit[i]=point1[2];

		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  vtxyfit[j]=point1[2];
		}else if(dist1>dist2){//point1 is correct
		  vtxxfit[i]=point2[0];
		  vtxzfit[i]=point2[1];
		  vtxyfit[i]=point2[2];

		  vtxxfit[j]=point2[0];
		  vtxzfit[j]=point2[1];
		  vtxyfit[j]=point2[2];
		}
		vtx_flag[i]=7;
		vtx_flag[j]=7;
	      }

	      mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;
	      mom_theta[j]=atan2(vtxzfit[j]-cz[j],vtxxfit[j]-cx[j])-acos(-1.)/2;

	      // vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      // vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      // vtxpxfit[j]=cos(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      // vtxpzfit[j]=sin(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      //	    G4cout<<"env_helm_field:"<<env_helm_field<<G4endl;
	    }
	  }
	  /////vertex reconstruction with beam
	  /*	  else if(sh_paID[i]==0. && sh_paID[i] != sh_paID[j]){
		  G4double rho1=rad[i];
		  G4double cx1=cx[i];

		  G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

		  double point1[2]={0};
		  G4int k;

		  if((cent_dist-(rho1+rho2))>0.){


		  G4double theta12=atan2(cz2-cz1,cx2-cx1);
		  G4double centr=rho1+(cent_dist-(rho1+rho2))/2;

		  point1[0]=cos(theta12)*centr+cx1;
		  point1[1]=sin(theta12)*centr+cz1;

		  G4double theta21=atan2(cz1-cz2,cx1-cx2);
		  G4double centr1=rho2+(cent_dist-(rho1+rho2))/2;
		  point2[0]=cos(theta21)*centr1+cx2;
		  point2[1]=sin(theta21)*centr1+cz2;

		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];
		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  }else  if((cent_dist+fmin(rho1,rho2))<fmax(rho1,rho2)){

		  if(rho1>=rho2){ //rho1>rho2
		  G4double theta12=atan2(cz2-cz1,cx2-cx1);
		  G4double centr=rho1-(rho1-cent_dist-rho2)/2; //rho1>rho2
		  point1[0]=cos(theta12)*centr+cx1;
		  point1[1]=sin(theta12)*centr+cz1;

		  G4double theta21=atan2(cz2-cz1,cx2-cx1);
		  G4double centr1=rho2+(rho1-cent_dist-rho2)/2.; //rho1>rho2
		  point2[0]=cos(theta21)*centr1+cx2;
		  point2[1]=sin(theta21)*centr1+cz2;
		  //		G4cout<<"test1"<<G4endl;

		  }else if(rho2>rho1){ //rho1<rho2
		  G4double theta12=atan2(cz1-cz2,cx1-cx2);
		  G4double centr=rho2-(rho2-cent_dist-rho1)/2; //rho1<rho2
		  point1[0]=cos(theta12)*centr+cx2;
		  point1[1]=sin(theta12)*centr+cz2;

		  G4double theta21=atan2(cz1-cz2,cx1-cx2);
		  G4double centr1=rho1+(rho2-cent_dist-rho1)/2; //rho1<rho2
		  point2[0]=cos(theta21)*centr1+cx1;
		  point2[1]=sin(theta21)*centr1+cz1;
		  }

		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];
		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  } else {

		  k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,point1,point2);
		  if(k == 0) {
		  G4cout << "no solution" << G4endl;
		  }else if(k>0){

		  G4double dist1=sqrt(pow(point1[0]-tpcData[i].tpcvtxx,2)+pow(point1[1]-tpcData[i].tpcvtxz,2));
		  G4double dist2=sqrt(pow(point2[0]-tpcData[i].tpcvtxx,2)+pow(point2[1]-tpcData[i].tpcvtxz,2));

		  if(dist1<=dist2){//point1 is correct
		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];

		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  }else if(dist1>dist2){//point1 is correct
		  vtxxfit[i]=point2[0];
		  vtxzfit[i]=point2[1];

		  vtxxfit[j]=point2[0];
		  vtxzfit[j]=point2[1];
		  }
		  }

		  mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;
		  mom_theta[j]=atan2(vtxzfit[j]-cz[j],vtxxfit[j]-cx[j])-acos(-1.)/2;

		  vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
		  vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
		  vtxpxfit[j]=cos(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
		  vtxpzfit[j]=sin(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);

		  }
		  }
	  */

	}
      }
    }


    /*
    //------------------&^-^-------------------//
    ////rungekutta study
    //------------------&^-^-------------------//
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //////////// rungekutta tracking based on LEPS TPC analyzer
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    int sector, lay;
    //    cir_r[kk]=rad[kk];
    //    cir_x[kk]=cx[kk];
    //    cir_z[kk]=cz[kk];
    double rkpar[5]={0};
    int iflag=0.;

    Switch sw;
    Track tracks[MAX_TRACK];
    initTrack(tracks);
    for(int kk=0; kk<tpctrNum; kk++){
    ///rkpara: x, y, u(px),v(py), q/p
    tracks[kk].rKInitPara[0]=x[kk][0];
    tracks[kk].rKInitPara[1]=y[kk][0];
    tracks[kk].rKInitPara[2]=vtxpxfit[kk]/sqrt(vtxpxfit[kk]*vtxpxfit[kk]+vtxpyfit[kk]*vtxpyfit[kk]+vtxpzfit[kk]*vtxpzfit[kk]);
    tracks[kk].rKInitPara[3]=vtxpyfit[kk]/sqrt(vtxpxfit[kk]*vtxpxfit[kk]+vtxpyfit[kk]*vtxpyfit[kk]+vtxpzfit[kk]*vtxpzfit[kk]);
    tracks[kk].rKInitPara[4]=tpcData[kk].tpcqq/sqrt(vtxpxfit[kk]*vtxpxfit[kk]+vtxpyfit[kk]*vtxpyfit[kk]+vtxpzfit[kk]*vtxpzfit[kk]);

    for(int j = 0; j < c[kk]; j++){
    tracks[kk].x[j][0] = x[kk][j];
    //	std::cout<<"j:"<<j<<"::"<<tracks[kk].x[j][0]<<std::endl;
    tracks[kk].x[j][1] = y[kk][j];
    tracks[kk].x[j][2] = z[kk][j];
    tracks[kk].numHits ++;
    //	std::cout<<"num_hits:"<<c[kk]<<std::endl;
    }
    }
    for(G4int kk=0; kk<tpctrNum; kk++){
    if(c[kk]>5.){
    RungeKuttaTracker rungekuttatrack(tracks+kk);
    }
    }

    */

    ///////////////////////vertex momentum for P_t
    G4int trn[MAX_TRACK];
    //// trancated mean --> now mean
    for(G4int i=0;i<tpctrNum;i++){
      trn[i]=c[i]*(truncated_mean_cut);
      G4double trtmp[MAX_TRACK]={0.000000000};
      for(G4int iii=0;iii<MAX_TRACK;iii++){
	trtmp[iii]=0.000000000;
      }

      for(G4int l=0;l<trn[i];l++){ //--> loop truncated number
	for(G4int k=0;k<c[i];k++){
	  if(l==0){
	    if(trtmp[l]<ede[i][k]){
	      trtmp[l]=ede[i][k];
	    }
	  }else if(l>0){
	    if(trtmp[l-1]>ede[i][k]){
	      if(trtmp[l]<ede[i][k]){
		trtmp[l]=ede[i][k];
	      }
	    }
	  }

	}//--loop end
      }
      for(G4int j=0;j<c[i];j++){
	if(trn[i]!=0.){
	  G4int sh_ch=1;
	  for(G4int jj=0;jj<trn[i];jj++){
	    if(ede[i][j]==trtmp[jj]){
	      sh_ch=-1;
	    }
	  }
	  if(sh_ch>0){
	    trmean[i]=trmean[i]+ede[i][j]/(c[i]-trn[i]);
	  }
	}else if(trn[i]==0.){
	  trmean[i]=trmean[i]+ede[i][j]/(c[i]-trn[i]);
	}

      }
      //      }
    }
    if(HitNum >= MaxHitsTPC){
      G4cerr << FUNC_NAME << " too much nhit (TPC) " << HitNum << G4endl;
    } else {

      for(G4int i=0; i<HitNum; i++){
	event.ntrk[event.nhittpc] = counterData[i].ntrk;
	hmap["Time"]->Fill(counterData[i].time);
	for(G4int j=0; j<G4ThreeVector::SIZE; ++j){
	  hmap[Form("Pos%d", j)]->Fill(counterData[i].pos[j]/CLHEP::mm);
	  hmap[Form("Mom%d", j)]->Fill(counterData[i].mom[j]/CLHEP::GeV);
	}
	event.xtpc[event.nhittpc] = counterData[i].pos[0]/CLHEP::mm;
	event.ytpc[event.nhittpc] = counterData[i].pos[1]/CLHEP::mm;
	event.ztpc[event.nhittpc] = counterData[i].pos[2]/CLHEP::mm;

	event.x0tpc[event.nhittpc] = counterData[i].pos0[0]/CLHEP::mm;
	event.y0tpc[event.nhittpc] = counterData[i].pos0[1]/CLHEP::mm;
	event.z0tpc[event.nhittpc] = counterData[i].pos0[2]/CLHEP::mm;

	event.resoX[event.nhittpc] = counterData[i].resoX;
	event.pxtpc[event.nhittpc] = counterData[i].mom[0]/CLHEP::GeV;
	event.pytpc[event.nhittpc] = counterData[i].mom[1]/CLHEP::GeV;
	event.pztpc[event.nhittpc] = counterData[i].mom[2]/CLHEP::GeV;
	event.pptpc[event.nhittpc] = sqrt(pow(counterData[i].mom[0], 2) +
                                          pow(counterData[i].mom[1], 2) +
                                          pow(counterData[i].mom[2], 2))/CLHEP::GeV;
	event.ititpc[event.nhittpc] = counterData[i].trackID;
	event.idtpc[event.nhittpc] = counterData[i].particleID;
	event.laytpc[event.nhittpc] = counterData[i].iLay;

	event.rowtpc[event.nhittpc] = counterData[i].iRow;
	event.iPadtpc[event.nhittpc] = padHelper::getPadID(event.laytpc[event.nhittpc], event.rowtpc[event.nhittpc]);
	TVector3 Point = padHelper::getPoint(event.iPadtpc[event.nhittpc]);
	event.xtpc_pad[event.nhittpc] = Point.x();
	event.ytpc_pad[event.nhittpc] = event.ytpc[event.nhittpc];
	event.ztpc_pad[event.nhittpc] = Point.z();

	event.dxtpc_pad[event.nhittpc] = event.x0tpc[event.nhittpc] - event.xtpc_pad[event.nhittpc];
	event.dytpc_pad[event.nhittpc] = event.y0tpc[event.nhittpc] - event.ytpc_pad[event.nhittpc];
	event.dztpc_pad[event.nhittpc] = event.z0tpc[event.nhittpc] - event.ztpc_pad[event.nhittpc];


	event.timetpc[event.nhittpc] = counterData[i].time/CLHEP::ns;
	event.betatpc[event.nhittpc] = counterData[i].beta;
	event.edeptpc[event.nhittpc] = counterData[i].edep/(CLHEP::MeV/CLHEP::mm);
	event.dedxtpc[event.nhittpc] = counterData[i].dedx;
	event.slengthtpc[event.nhittpc] = counterData[i].slength/CLHEP::mm;
	event.tlengthtpc[event.nhittpc] = counterData[i].tlength/CLHEP::mm;
	event.nthlay[event.nhittpc] = counterData[i].iLay;
	event.nthpad[event.nhittpc] = counterData[i].iPad;
	event.laypad[event.nhittpc][event.nthlay[event.nhittpc]][event.nthpad[event.nhittpc]]
	  = event.laypad[event.nhittpc][event.nthlay[event.nhittpc]][event.nthpad[event.nhittpc]]+1.;
	event.nhittpc += 1;
      }
    }
    //
    // TPC
    //
    for(G4int i=0; i<tpctrNum; i++){
      //      G4cout<<"abs"<<abs(env_helm_field)<<G4endl;
      //      G4cout<<"fabs"<<fabs(env_helm_field)<<G4endl;
      // anaRoot.FillTPCData(tpcData[i].tpcpx,
      // 			  tpcData[i].tpcpy,tpcData[i].tpcpz,
      // 			  tpcData[i].tpcpp,
      // 			  tpcData[i].tpcpid, tpcData[i].tpcparentid, tpcData[i].tpcparentid_pid,
      // 			  tpcData[i].tpcqq,
      // 			  tpcData[i].tpcpm,tpcData[i].tpcde,
      // 			  tpcData[i].tpclen,mean[i],trmean[i],
      // 			  tpcData[i].tpclay,
      // 			  tpcData[i].tpcvtxpx,tpcData[i].tpcvtxpy,tpcData[i].tpcvtxpz,
      // 			  tpcData[i].tpcvtxx,tpcData[i].tpcvtxy,tpcData[i].tpcvtxz,
      // 			  vtxxfit[i],vtxyfit[i],vtxzfit[i],
      // 			  vtxpxfit[i],Pz[i],vtxpzfit[i],cir_r[i]*(0.299792458)*fabs(env_helm_field),
      // 			  cir_r[i],cir_x[i],cir_z[i],test[i],
      // 			  vtx_flag[i], a_fory[i], b_fory[i]
      // 			 );
    }
  }//trigger parts

  m_tree->Fill();

  return 0;
}

//_____________________________________________________________________________
void
AnaManager::SetHitData(const VHitInfo* hit)
{
  if(hit && hit->GetParticle()){
    const auto& name = hit->GetDetectorName();
    event.hits[name].push_back(*hit->GetParticle());
  }
}

//_____________________________________________________________________________
void
AnaManager::SetCounterData(G4int ntrk, G4double time, G4ThreeVector pos,
                           G4ThreeVector mom,
                           G4int track, G4int particle,
                           G4int iLay,  G4int iRow, G4double beta,
                           G4double edep, G4int parentid,
                           G4double tlength, G4double slength)
{
  G4int hitnum = HitNum;
  G4bool flag=true;
  if (hitnum >= MaxTrack) {
    fprintf(stderr, "AnaManager::SetCounterData Too Much multiplicity %d\n",
            hitnum);
    return;
  }

  //  G4ThreeVector tar_pos(0.,0.*CLHEP::mm,-150.*CLHEP::mm);

  G4ThreeVector tar_pos(0.,0., -143);
  G4ThreeVector sh_pos(0.,0.,0.);
  sh_pos=pos-tar_pos;

  G4double sh_r = sh_pos.r();
  G4double sh_theta = sh_pos.theta();
  G4double sh_phi = sh_pos.phi();

  G4double sh_x = sh_r*sin(sh_theta)*cos(sh_phi);
  G4double sh_y = sh_r*sin(sh_theta)*sin(sh_phi);
  G4double sh_z = sh_r*cos(sh_theta);

  counterData[hitnum].particleID = particle;
  ////shhwang check, check a multiplicity of layers
  for(G4int i=0;i<hitnum;i++){
    if((counterData[i].iLay == iLay &&
        counterData[i].ntrk == ntrk)){
      flag = false;
    }
  }
  flag=true;
  if(flag == true){
    counterData[hitnum].ntrk = ntrk;
    counterData[hitnum].time = time;
    counterData[hitnum].beta = beta;
    counterData[hitnum].dedx = edep/slength;
    counterData[hitnum].edep = edep;
    counterData[hitnum].slength = slength;
    counterData[hitnum].tlength = tlength;

    //////shhwang position smearing////


    G4double sh_alpha =  atan2(sh_x,sh_z);
    G4double sh_rho =  sqrt(pow(sh_z,2)+pow(sh_x,2));
    //    G4double sh_dalpha = 0.300*CLHEP::mm/sh_rho; // rho * theta = arc --> sigma=300 um
    //    G4double sh_smear_alpha = CLHEP::RandGauss::shoot(sh_alpha, sh_dalpha);
    G4double sh_sigmaY = 0.500*CLHEP::mm; //--> smearing : 400 um

    G4double ang_sh=atan2(sh_pos.getY(),sh_pos.getX());
    // G4double ang_check=0;
    if(ang_sh>acos(-1.)){
      ang_sh=ang_sh-2*acos(-1.);
      // ang_check=-1;
    }else{
      // ang_check=1.;
    }
    // G4double arc_sh=pad_in[iLay]*ang_sh;
    // G4int ith_pad_in=arc_sh/pad_in_width;
    // G4int ith_pad_out=arc_sh/pad_out_width;


    // G4double delta_x=0.;

    // if(m_pad_config ==1){
    //   if(iLay<pad_in_num){   // const G4int pad_in_num=10;
    // 	delta_x=arc_sh-ith_pad_in*pad_in_width;
    // 	if(delta_x<0){
    // 	  delta_x=delta_x+pad_in_width;
    // 	}
    //   }else if(iLay>=pad_in_num){
    // 	delta_x=arc_sh-ith_pad_out*pad_out_width;
    // 	if(delta_x<0){
    // 	  delta_x=delta_x+pad_out_width;
    // 	}
    //   }
    // }else if(m_pad_config ==2){
    //   if(iLay<pad_in_num){   //

    // 	ith_pad_in=arc_sh/seg_width[iLay];
    // 	delta_x=arc_sh-ith_pad_in*seg_width[iLay];
    // 	if(delta_x<0){
    // 	  delta_x=delta_x+seg_width[iLay];
    // 	}

    //   }else if(iLay>=pad_in_num){
    // 	ith_pad_out=arc_sh/seg_width[iLay];
    // 	delta_x=arc_sh-ith_pad_out*seg_width[iLay];
    // 	if(delta_x<0){
    // 	  delta_x=delta_x+seg_width[iLay];
    // 	}

    //   }
    // }


    ///include saho-san's code

    // G4double compy=0.;
    G4double compx=0.;

    // G4int n_electron=0; //output
    // G4int n_pad=0;//output
    // G4double x_rms=0;//output
    // G4double x_track = delta_x;// pad edge
    // G4double y_track = double(250.+50.+pos.getY());// 50 is office set
    // if(y_track<0 || y_track>600) G4cout<<"Y_track estimation is wrong(y_track:y_pos)"<<y_track<<":"<<pos.getY()<<G4endl;


    // G4double dxdz_track=tan(acos((sh_pos.getX()*mom.getX()+sh_pos.getZ()*mom.getZ())/(sqrt(sh_pos.getX()*sh_pos.getX()+sh_pos.getZ()*sh_pos.getZ())*sqrt(mom.getX()*mom.getX()+mom.getZ()*mom.getZ()))));


    // //pad and particle angle


    // G4double dydz_track=0;

    // if(m_pad_config ==1){
    //   if(iLay<pad_in_num){
    // 	ResHypTPC reshyptpc(pad_in_width, pad_length_in, 0.1,0.18, 0);
    // 	compx = reshyptpc.getXDeviation(n_electron, n_pad, x_rms, x_track, y_track, dxdz_track, dydz_track);
    // 	compy = reshyptpc.getYDeviation(pos.getY());
    //   }else if(iLay>=pad_in_num){
    // 	ResHypTPC reshyptpc(pad_out_width, pad_length_out, 0.1,0.18, 0);
    // 	compx = reshyptpc.getXDeviation(n_electron, n_pad, x_rms, x_track, y_track, dxdz_track, dydz_track);
    // 	compy = reshyptpc.getYDeviation(pos.getY());
    //   }


    // }else if(m_pad_config ==2){
    //   if(iLay<pad_in_num){
    // 	ResHypTPC reshyptpc(double(seg_width[iLay]-0.5), pad_length_in, 0.1,0.18, 0);
    // 	compx = reshyptpc.getXDeviation(n_electron, n_pad, x_rms, x_track, y_track, dxdz_track, dydz_track);
    // 	compy = reshyptpc.getYDeviation(pos.getY());
    //   }else if(iLay>=pad_in_num){
    //ResHypTPC reshyptpc(double(seg_width[iLay]-0.5), pad_length_out, 0.1,0.18, 0);
    // 	compx = reshyptpc.getXDeviation(n_electron, n_pad, x_rms, x_track, y_track, dxdz_track, dydz_track);
    // 	compy = reshyptpc.getYDeviation(pos.getY());
    //   }
    // }
    compx = GetTransverseRes(sh_y);
    double s_compx = CLHEP::RandGauss::shoot(0.,compx);

    // std::cout<<"compx="<<compx<<", sh_sigmaY"<<sh_sigmaY<<std::endl;
    // getchar();
    //G4double sh_dalpha = compx/sh_rho; // rho * theta = arc --> from sako-san's code
    //      G4double sh_dalpha = s_compx/sh_rho; // rho * theta = arc --> from sako-san's code
    G4double sh_dalpha = atan2(s_compx, sh_rho); // rho * theta = arc --> from sako-san's code
    G4double sh_smear_alpha = sh_alpha+sh_dalpha;
    //    G4cout<<compx<<":"<<sh_dalpha<<G4endl;
    // G4double randx = sh_rho*(sin(sh_smear_alpha)-sin(sh_alpha));
    // G4double randz = sh_rho*(cos(sh_smear_alpha)-cos(sh_alpha));
    //G4double s0 = 0.204;// mm HIMAC result
    // G4double s0 = 0.199;// mm HIMAC result
    // G4double randx_com = CLHEP::RandGauss::shoot(0.,s0);
    // G4double randz_com = CLHEP::RandGauss::shoot(0.,s0);
    // G4double smear_x = sqrt(randx*randx + randx_com*randx_com)*(randx_com/fabs(randx_com));
    // G4double smear_z = sqrt(randz*randz + randz_com*randz_com)*(randz_com/fabs(randz_com));

    counterData[hitnum].resoX = compx;

    // counterData[hitnum].pos[G4ThreeVector::Z] = sh_rho*cos(sh_alpha)+smear_z+tar_pos.getZ();
    // counterData[hitnum].pos[G4ThreeVector::X] = sh_rho*sin(sh_alpha)+smear_x;
    counterData[hitnum].pos[G4ThreeVector::Z] = sh_rho*cos(sh_smear_alpha)+tar_pos.getZ();
    counterData[hitnum].pos[G4ThreeVector::X] = sh_rho*sin(sh_smear_alpha);
    counterData[hitnum].pos[G4ThreeVector::Y] = CLHEP::RandGauss::shoot(sh_y,sh_sigmaY);

    counterData[hitnum].pos0[G4ThreeVector::X] = pos.getX();
    counterData[hitnum].pos0[G4ThreeVector::Y] = pos.getY();
    counterData[hitnum].pos0[G4ThreeVector::Z] = pos.getZ();

    counterData[hitnum].mom[G4ThreeVector::X] = mom.getX();
    counterData[hitnum].mom[G4ThreeVector::Y] = mom.getY();
    counterData[hitnum].mom[G4ThreeVector::Z] = mom.getZ();

    counterData[hitnum].trackID = track;
    counterData[hitnum].particleID = particle;
    counterData[hitnum].iLay = iLay;
    G4int iPad=0.;

    if(m_pad_config == 2){
      //    G4bool pass_check=false;
      G4bool pass_check=true;
      G4double cur_angle= (acos(-1.)-atan2(sh_x,sh_z))*180./acos(-1.);
      //    G4cout<<"--------------------"<<G4endl;
      //    G4cout<<"currrent angle:"<<cur_angle<<G4endl;
      //    G4cout<<"layer angle:"<<angle[iLay]<<G4endl;
      //    G4cout<<"seg angle:"<<seg_angle[iLay]<<G4endl;
      if(iLay<pad_in_num){
        G4double check_num_pads=(cur_angle)/seg_angle[iLay];
        //      G4cout<<check_num_pads<<G4endl;
        iPad=int(check_num_pads);
      }else if(iLay>=pad_in_num){
        G4double check_num_pads=(cur_angle-angle[iLay])/seg_angle[iLay];
        //      G4cout<<check_num_pads<<G4endl;
        iPad=int(check_num_pads);
      }
      if(iPad>numpads[iLay]){
        G4cout<<"this code has a error(iPad:numpads)-->"<<iPad<<":"<<numpads[iLay]<<G4endl;
      }
      if(pass_check){
        //    G4cout<<"iLay:"<<iLay<<G4endl;
        //        G4cout<<"iLay:"<<iLay<<G4endl;
        //        G4cout<<"iPad:"<<iPad<<G4endl;
        counterData[hitnum].iPad = iPad;
      }else{
        G4cout<<"wrong:"<<iLay<<G4endl;
      }
    }

    counterData[hitnum].iRow = iRow;
    counterData[hitnum].parentID = parentid;
    HitNum++;

    if(particle==321)
      HitNum_K++;

    if(particle==2212)
      HitNum_p++;
  }

  return;
}

//_____________________________________________________________________________
void
AnaManager::SetFermiMomentum(const G4ThreeVector& p)
{
  for(G4int i=0; i<G4ThreeVector::SIZE; ++i){
    TString key = Form("Fermi%d", i);
    hmap[key]->Fill(p[i]*CLHEP::GeV);
  }
}

//_____________________________________________________________________________
void
AnaManager::SetTPCData(G4int tpctr2, G4int tpcpid2, G4int tpcparentid2,
                       G4int tpcparentid_pid2, G4double tpcpx2,
                       G4double tpcpy2, G4double tpcpz2,
                       G4double /* tpcpp2 */,
                       G4int tpcqq2, G4double tpcpm2, G4double tpcde2,
                       G4double tpclen2, G4int tpclay2,
                       G4double vtxpxtpc2, G4double vtxpytpc2,
                       G4double vtxpztpc2,
                       G4double vtxxtpc2, G4double vtxytpc2,
                       G4double vtxztpc2, G4double vtxene2)
{
  G4int hitnum = tpctr2;

  // G4double theta=acos(tpcpz2/tpcpp2);

  tpcData[hitnum].tpctr = tpctr2;
  tpcData[hitnum].tpcpid = tpcpid2;
  tpcData[hitnum].tpcparentid = tpcparentid2;
  tpcData[hitnum].tpcparentid_pid = tpcparentid_pid2;

  //// w/o smearing
  tpcData[hitnum].tpcpx = tpcpx2;
  tpcData[hitnum].tpcpy = tpcpy2;
  tpcData[hitnum].tpcpz = tpcpz2;
  //// with smearing
  //    tpcData[hitnum].tpcpx = px;
  //    tpcData[hitnum].tpcpz = pz;
  //    tpcData[hitnum].tpcpy = py;


  //kine E = sqrt(p^2+m^2)-m
  //p=sqrt((E+m)^2-m^2)
  G4double totalmom=sqrt(pow(vtxene2+tpcpm2,2)-pow(tpcpm2,2));
  tpcData[hitnum].tpcvtxpx = totalmom*vtxpxtpc2;
  tpcData[hitnum].tpcvtxpy = totalmom*vtxpytpc2;
  tpcData[hitnum].tpcvtxpz = totalmom*vtxpztpc2;

  tpcData[hitnum].tpcvtxx = vtxxtpc2;
  tpcData[hitnum].tpcvtxy = vtxytpc2;
  tpcData[hitnum].tpcvtxz = vtxztpc2;

  //// with smearing
  //  tpcData[hitnum].tpcpp = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

  //// w/o smearing
  tpcData[hitnum].tpcpp = sqrt(pow(tpcpx2,2)+pow(tpcpy2,2)+pow(tpcpz2,2));
  //  tpcData[hitnum].tpcppfit = sqrt(pow(,2)+pow(tpcpy2,2));

  tpcData[hitnum].tpcqq = tpcqq2;
  tpcData[hitnum].tpcpm = tpcpm2;
  tpcData[hitnum].tpclen = tpclen2;
  tpcData[hitnum].tpcdedx = tpcde2;
  tpcData[hitnum].tpclay = tpclay2;
  tpctrNum++;
  return;
}

//_____________________________________________________________________________
void
AnaManager::SetPrimaryParticle(G4int id, G4int pdg,
                               const G4LorentzVector& p,
                               const G4LorentzVector& v,
                               G4bool is_virtual_beam)
{
  G4int id1 = is_virtual_beam ? -1 : 1;
  G4int id2 = id;
  for (const auto& ptcl: event.hits["PRM"]) {
    if (ptcl.GetMother(0) == id1 && ptcl.GetMother(1) == id2) {
      G4cerr << FUNC_NAME << " id1=" << id1 << ", id2=" << id2
             << " is already set" << G4endl;
    }
  }
  TParticle particle(pdg,
                     0, // fStatus
                     id1, // fMother[0]
                     id2, // fMother[1]
                     0, // fDaughter[0]
                     0, // fDaughter[1]
                     TLorentzVector(p.px(), p.py(), p.pz(), p.e()),
                     TLorentzVector(v.x(), v.y(), v.z(), v.t()));
  event.hits["PRM"].push_back(particle);
}

//_____________________________________________________________________________
void
AnaManager::SetGeneratorID(G4int generator)
{
  event.generator = generator;
}

//_____________________________________________________________________________
void
AnaManager::SetModeID(G4int mode)
{
  event.mode = mode;
}

//_____________________________________________________________________________
void
AnaManager::SetIncID(G4int inc)
{
  event.inc = inc;
}

/*************************************
 *************************************/
void initTrack(Track* tracks){
  static const std::string funcname = "[InitTrack]";
  int i,j;
  //  G4cout<<"init track"<<G4endl;
  for(i = 0; i < MAX_TRACK; i++){
    tracks[i].nout   =  0;
    tracks[i].ngood  =  0;
    tracks[i].igroup = -1;
    tracks[i].trkQual = -1;
    tracks[i].numHits = 0;
    ////    tracks[i].numSectors = 0;
    tracks[i].numLayers = 0;
    tracks[i].charge = 1000;
    tracks[i].totalLength = 0.0;
    tracks[i].totalLengthTOF = 0.0; /*NTPC TOF*/
    tracks[i].meanAdc = -1.0;
    tracks[i].chi2Pad = -1.0;
    tracks[i].chi2Z       = -1.0;
    tracks[i].chi2Prob = -1.0;
    tracks[i].chi2 = -1.0;
    tracks[i].radius = 0.0;
    tracks[i].center[0] = tracks[i].center[1] = 1000;
    tracks[i].rKNumIter = -1;
    tracks[i].CrossOuter = -1;
    tracks[i].mom[0] = tracks[i].mom[1] =
      tracks[i].mom[2] = 10000.0;
    tracks[i].mom[3] = -1.0;
    tracks[i].resVirtual[0] = tracks[i].resVirtual[1]
      =tracks[i].resVirtual[2] = tracks[i].resVirtual[3] = -100;
    tracks[i].xOnTrack[0] = tracks[i].xOnTrack[1]
      = tracks[i].xOnTrack[2] = 1000.0;
    tracks[i].RKPFinal[0] =
      tracks[i].RKPFinal[1] = tracks[i].RKPFinal[2]  = -1000;
    //    for(j=0 ; j<NUM_SECTOR ; j++){
    //      tracks[i].numLayersinSec[j] = 0;
    //    }
    for(j=0; j < MAX_ITERATION; j++){
      tracks[i].rKChi2[j] = 1000.;
    }

    for(j=0; j < NUM_PARA_RK; j++){
      tracks[i].rKInitPara[j] =
        tracks[i].rKFinalPara[j] = -1000.;
    }

    /*  tracks[i].hitPattern = 0;*/
    for(j = 0; j < MAX_HIT_IN_TRACK; j++){
      tracks[i].ibad[j]  = 10;
      tracks[i].zbad[j]  = 10;
      //      tracks[i].hit[j]    = NULL;
      tracks[i].resPad[j] = 1000;
      tracks[i].resPady[j] = 1000; /*NTPC*/
      tracks[i].resZ[j]   = 1000;
      tracks[i].sector[j] = 100;
      tracks[i].lay[j] = 100;
      tracks[i].x[j][0] = 1000.;/*originally 100.*/
      tracks[i].x[j][1] = 1000.;/*originally 100.*/
      tracks[i].x[j][2] = -10000.;
      tracks[i].err[j][0] = tracks[i].err[j][1]
        = tracks[i].err[j][2] = 1000.;
      tracks[i].arcLen[j] = -1.0;
      tracks[i].resPad[j] = tracks[i].resPady[j] = tracks[i].resZ[j] =
        tracks[i].rKresXYZ[j][0] = tracks[i].rKresXYZ[j][1] =
        tracks[i].rKresXYZ[j][2] = tracks[i].rKresXYZ[j][3] =
        tracks[i].initRes[j][0] = tracks[i].initRes[j][1] =
        tracks[i].initRes[j][2] = tracks[i].phi_local[j]
        = -10000.;

    }
  }
}


void initTrack_ku(Track* tracks){
  static const std::string funcname = "[InitTrack]";
  int i,j;
  //  G4cout<<"init track"<<G4endl;
  for(i = 0; i < MAX_TRACK; i++){
    tracks[i].nout   =  0;
    tracks[i].ngood  =  0;
    tracks[i].igroup = -1;
    tracks[i].trkQual = -1;
    tracks[i].numHits = 0;
    ////    tracks[i].numSectors = 0;
    tracks[i].numLayers = 0;
    tracks[i].charge = 1000;
    tracks[i].totalLength = 0.0;
    tracks[i].totalLengthTOF = 0.0; /*NTPC TOF*/
    tracks[i].meanAdc = -1.0;
    tracks[i].chi2Pad = -1.0;
    tracks[i].chi2Z       = -1.0;
    tracks[i].chi2Prob = -1.0;
    tracks[i].chi2 = -1.0;
    tracks[i].radius = 0.0;
    tracks[i].center[0] = tracks[i].center[1] = 1000;
    tracks[i].rKNumIter = -1;
    tracks[i].CrossOuter = -1;
    tracks[i].mom[0] = tracks[i].mom[1] =
      tracks[i].mom[2] = 10000.0;
    tracks[i].mom[3] = -1.0;
    tracks[i].resVirtual[0] = tracks[i].resVirtual[1]
      =tracks[i].resVirtual[2] = tracks[i].resVirtual[3] = -100;
    tracks[i].xOnTrack[0] = tracks[i].xOnTrack[1]
      = tracks[i].xOnTrack[2] = 1000.0;
    tracks[i].RKPFinal[0] =
      tracks[i].RKPFinal[1] = tracks[i].RKPFinal[2]  = -1000;
    //    for(j=0 ; j<NUM_SECTOR ; j++){
    //      tracks[i].numLayersinSec[j] = 0;
    //    }
    for(j=0; j < MAX_ITERATION; j++){
      tracks[i].rKChi2[j] = 1000.;
    }

    for(j=0; j < NUM_PARA_RK; j++){
      tracks[i].rKInitPara[j] =
        tracks[i].rKFinalPara[j] = -1000.;
    }

    /*  tracks[i].hitPattern = 0;*/
    for(j = 0; j < MAX_HIT_IN_TRACK; j++){
      tracks[i].ibad[j]  = 10;
      tracks[i].zbad[j]  = 10;
      //      tracks[i].hit[j]    = NULL;
      tracks[i].resPad[j] = 1000;
      tracks[i].resPady[j] = 1000; /*NTPC*/
      tracks[i].resZ[j]   = 1000;
      tracks[i].sector[j] = 100;
      tracks[i].lay[j] = 100;
      tracks[i].x[j][0] = 1000.;/*originally 100.*/
      tracks[i].x[j][1] = 1000.;/*originally 100.*/
      tracks[i].x[j][2] = -10000.;
      tracks[i].res[j] = 0.;
      tracks[i].err[j][0] = tracks[i].err[j][1]
        = tracks[i].err[j][2] = 1000.;
      tracks[i].arcLen[j] = -1.0;
      tracks[i].resPad[j] = tracks[i].resPady[j] = tracks[i].resZ[j] =
        tracks[i].rKresXYZ[j][0] = tracks[i].rKresXYZ[j][1] =
        tracks[i].rKresXYZ[j][2] = tracks[i].rKresXYZ[j][3] =
        tracks[i].initRes[j][0] = tracks[i].initRes[j][1] =
        tracks[i].initRes[j][2] = tracks[i].phi_local[j]
        = -10000.;

    }
  }
}

/*************************************
 *************************************/
/*void setTrack(Track* tracks,int ntrk, double* x, double* y, double* z, double* ede, double* nhit){
  static const std::string funcname = "[SetTrack]";
  int i,j;
  for(i = 0; i < ntrk; i++){
  for(j = 0; j < nhit[i]; j++){
  tracks[i].lay[j] = 100;
  tracks[i].x[j][0] = x[j];
  tracks[i].x[j][1] = y[j];
  tracks[i].x[j][2] = z[j];
  }
  }
  }
*/
/*************************************
 *************************************/
