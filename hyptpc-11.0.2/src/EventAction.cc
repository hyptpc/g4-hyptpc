// -*- C++ -*-

#include "EventAction.hh"

#include <fstream>

#include <G4RunManager.hh>
#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4VVisManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4VTrajectory.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>

#include "AnaManager.hh"
#include "TPCPadHelper.hh"
#include "BACSD.hh"
#include "ConfMan.hh"
#include "DetectorConstruction.hh"
#include "FuncName.hh"
#include "BH2SD.hh"
#include "FTOFSD.hh"
#include "HTOFSD.hh"
#include "BACSD.hh"
#include "TPCSD.hh"
#include "TPCPadSD.hh"
#include "SCHSD.hh"
#include "SDCSD.hh"
#include "TargetSD.hh"
#include "VPSD.hh"
#include "KVCSD.hh"

#define SDCLASS(x) x ## Hit

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
EventAction::EventAction()
  : G4UserEventAction()
{
}

//_____________________________________________________________________________
EventAction::~EventAction()
{
}

//_____________________________________________________________________________
void
EventAction::BeginOfEventAction(const G4Event*)
{
#ifdef DEBUG
  G4cout << FUNC_NAME << G4endl;
#endif
  gAnaMan.BeginOfEventAction();
}

//_____________________________________________________________________________
void
EventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4int eventID = anEvent->GetEventID();
  if (eventID % 100 == 0) {
    G4cout << FUNC_NAME << G4endl
	   << "   Event number = " << eventID << G4endl;
  }
  auto SDManager = G4SDManager::GetSDMpointer();
  auto HCTE = anEvent->GetHCofThisEvent();
  if (!HCTE) return;

  {
    static const auto id = SDManager->GetCollectionID("BH2/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<BH2Hit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("BH2", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("TGT/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<TargetHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("TGT", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("HTOF/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<HTOFHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("HTOF", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("FTOF/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<FTOFHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("FTOF", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("BAC/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<BACHit>*>(HCTE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("BAC", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("KVC/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<KVCHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("KVC", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("VP/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<VPHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("VP", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("TPC/hit");
    if(id > 0){
      auto HC = dynamic_cast<G4THitsCollection<TPCHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("TPC", HC->entries());
    }
  }

  {
    static const G4int id = SDManager->GetCollectionID("TPCPad/hit");
    if( id > 0 ){
      auto HC = (G4THitsCollection<TPCPadHit>*)(HCTE->GetHC(id));
      G4int nhits= HC -> entries();
      G4int pidtr[MaxHitsTPC]={0};
      G4int ptidtpc[MaxHitsTPC]={0};
      G4int ptidtpc_pid[MaxHitsTPC]={0};
      G4double pmtpc[MaxHitsTPC]={0};
      G4int qqtpc[MaxHitsTPC]={0};
      G4double pxtpc[MaxHitsTPC]={0};
      G4double pytpc[MaxHitsTPC]={0};
      G4double pztpc[MaxHitsTPC]={0};
      G4double pptpc[MaxHitsTPC]={0};
      G4double vtxxtpc[MaxHitsTPC]={0};
      G4double vtxytpc[MaxHitsTPC]={0};
      G4double vtxztpc[MaxHitsTPC]={0};
      G4double vtxpxtpc[MaxHitsTPC]={0};
      G4double vtxpytpc[MaxHitsTPC]={0};
      G4double vtxpztpc[MaxHitsTPC]={0};
      // G4double vtxpptpc[MaxHitsTPC]={0};
      G4double vtxenetpc[MaxHitsTPC]={0};
      G4double detpc[MaxHitsTPC]={0};
      G4int laytpc[MaxHitsTPC]={0};
      G4double lentpc[MaxHitsTPC]={0};
      G4int nparticle=0;
      // G4cout << "TPC  " << nhits << G4endl;
      for( G4int i=0; i<nhits; ++i ){
	if(nparticle >20) continue;
	if(TPCPadHelper::GetDeadCon((*HC)[i]->GetPadLay(), (*HC)[i]->GetPadRow())){
	  continue;
	  std::cout<<"dead pad!!"<<std::endl;
	}

	G4ThreeVector vtxpos = (*HC)[i]-> GetVtxPosition();
	G4ThreeVector vtxmom = (*HC)[i]-> GetVtxMomentum();
	G4double vtxene =(*HC)[i]-> GetVtxEnergy();
	G4ThreeVector xyz = (*HC)[i]-> GetPosition();
	G4ThreeVector mom = (*HC)[i]-> GetMomentum();
	G4double tof= (*HC)[i]-> GetTOF();
	G4int tid = (*HC)[i]-> GetTrackID();
	G4int ptid = (*HC)[i]-> GetParentID();
	G4int ptid_pid = (*HC)[i]-> GetParentID_pid();
	G4int pid = (*HC)[i]-> GetParticleID();
	G4double mass = (*HC)[i]-> GetMass();
	G4int charge = (*HC)[i]-> GetCharge();
	// std::cout<<"pid="<<pid<<", mass="<<mass<<", charge="<<charge<<std::endl;
	// getchar();
	G4int ilay = (*HC)[i]-> GetPadLay();
	//      G4double mass = (*HC)[i]-> GetPDGMass(); //mass(GeV)
	G4int parentid = (*HC)[i]-> GetParentID();
	G4int parentpid = (*HC)[i]-> GetParentID_pid();
	G4double tlength = (*HC)[i]-> GettLength();
	G4int irow=(*HC)[i]-> GetPadRow();
	G4double beta = (*HC)[i]-> GetBeta();
	G4double edep = (*HC)[i]-> GetEdep();
	G4double slength = (*HC)[i]-> GetsLength();
	//      G4VTrajectoryPoint *tp -> HC->GetPoint(i);
	//    G4int nhits= HC -> entries();
	if( nparticle==0 ){
	  qqtpc[nparticle]=charge;
	  pmtpc[nparticle]=mass;
	  detpc[nparticle]=detpc[nparticle]+edep;
	  // vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));
	  if(ilay>-1){
	    laytpc[nparticle]=laytpc[nparticle]+1;
	  }
	  pidtr[nparticle]=pid;
	  pxtpc[nparticle]=mom[0];
	  pytpc[nparticle]=mom[1];
	  pztpc[nparticle]=mom[2];

	  //////////////////////vertex information /////////////////////////
	  vtxpxtpc[nparticle]=vtxmom[0];
	  vtxpytpc[nparticle]=vtxmom[1];
	  vtxpztpc[nparticle]=vtxmom[2];
	  // vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

	  vtxxtpc[nparticle]=vtxpos[0];
	  vtxytpc[nparticle]=vtxpos[1];
	  vtxztpc[nparticle]=vtxpos[2];

	  vtxenetpc[nparticle]=vtxene;

	  pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	  lentpc[nparticle]=tlength;
	  ptidtpc[nparticle]=ptid;
	  ptidtpc_pid[nparticle]=ptid_pid;
	  nparticle=nparticle+1;

	}else if( nparticle>0 ){
	  //	G4cout<<nparticle<<G4endl;
	}
	if( (pidtr[nparticle-1] != pid) || (pidtr[nparticle-1] == pid && vtxpxtpc[nparticle-1] != vtxmom[0] && vtxpytpc[nparticle-1] != vtxmom[1] && vtxpztpc[nparticle-1] != vtxmom[2])){
	  qqtpc[nparticle]=charge;
	  pmtpc[nparticle]=mass;
	  detpc[nparticle]=detpc[nparticle]+edep;
	  if(ilay>-1){
	    laytpc[nparticle]=laytpc[nparticle]+1;
	  }
	  pidtr[nparticle]=pid;
	  pxtpc[nparticle]=mom[0];
	  pytpc[nparticle]=mom[1];
	  pztpc[nparticle]=mom[2];

	  vtxpxtpc[nparticle]=vtxmom[0];
	  vtxpytpc[nparticle]=vtxmom[1];
	  vtxpztpc[nparticle]=vtxmom[2];
	  // vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));
	  vtxenetpc[nparticle]=vtxene;

	  vtxxtpc[nparticle]=vtxpos[0];
	  vtxytpc[nparticle]=vtxpos[1];
	  vtxztpc[nparticle]=vtxpos[2];

	  pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	  lentpc[nparticle]=tlength;
	  ptidtpc[nparticle]=ptid;
	  ptidtpc_pid[nparticle]=ptid_pid;
	  nparticle=nparticle+1;
	}else if (pidtr[nparticle-1] == pid && vtxpxtpc[nparticle-1] == vtxmom[0] && vtxpytpc[nparticle-1] == vtxmom[1] && vtxpztpc[nparticle-1] == vtxmom[2]){
	  if( ptidtpc[nparticle-1] != ptid){
	    qqtpc[nparticle]=charge;
	    pmtpc[nparticle]=mass;
	    detpc[nparticle]=detpc[nparticle]+edep;
	    if(ilay>-1){
	      laytpc[nparticle]=laytpc[nparticle]+1;
	    }
	    pidtr[nparticle]=pid;
	    pxtpc[nparticle]=mom[0];
	    pytpc[nparticle]=mom[1];
	    pztpc[nparticle]=mom[2];

	    vtxpxtpc[nparticle]=vtxmom[0];
	    vtxpytpc[nparticle]=vtxmom[1];
	    vtxpztpc[nparticle]=vtxmom[2];
	    // vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

	    vtxenetpc[nparticle]=vtxene;

	    vtxxtpc[nparticle]=vtxpos[0];
	    vtxytpc[nparticle]=vtxpos[1];
	    vtxztpc[nparticle]=vtxpos[2];

	    pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	    lentpc[nparticle]=tlength;
	    ptidtpc[nparticle]=ptid;
	    ptidtpc_pid[nparticle]=ptid_pid;
	    nparticle=nparticle+1;
	  }
	  else{
	    detpc[nparticle-1]=detpc[nparticle-1]+edep;
	    if(ilay>-1){
	      laytpc[nparticle-1]=laytpc[nparticle-1]+1;
	    }
	    lentpc[nparticle-1]=tlength;
	  }
	}
	//      if(ilay>-1){ //--> ilay 1 : target ilay 0 : TPC, layer is from 2 to 38.
	if(ilay>-1){ //-->  -1 : TPC, layer is from 0 to 38. 2012.10.30
	  int restype = gConf.Get<G4int>("ResType");
	  switch(restype){
	  case 0:
	    gAnaMan.SetCounterDataSimple( nparticle-1,tof, xyz, mom, tid, pid, ilay,
					  irow, beta, edep/CLHEP::MeV, parentid, tlength,slength );
	    break;
	  case 1:
	    gAnaMan.SetCounterDataExp( nparticle-1,tof, xyz, mom, tid, pid, ilay,
				       irow, beta, edep/CLHEP::MeV, parentid, tlength,slength );
	    break;
	  default:
	    gAnaMan.SetCounterDataSimple( nparticle-1,tof, xyz, mom, tid, pid, ilay,
					  irow, beta, edep/CLHEP::MeV, parentid, tlength,slength );
	    G4cout<<"TPC Resolution type is not determined, use constant resolution"<<G4endl;
	    break;
	  }
	
	}
      }
      for(G4int i=0;i<nparticle;i++){
	gAnaMan.SetTPCData(i, pidtr[i], ptidtpc[i], ptidtpc_pid[i], pxtpc[i],pytpc[i],pztpc[i],pptpc[i], qqtpc[i], pmtpc[i], detpc[i], lentpc[i],laytpc[i],
			   vtxpxtpc[i],vtxpytpc[i],vtxpztpc[i],
			   vtxxtpc[i],vtxytpc[i],vtxztpc[i], vtxenetpc[i]);
      }
    }
  }

#if 0
  auto trajectoryContainer = anEvent->GetTrajectoryContainer();
  if(trajectoryContainer && G4VVisManager::GetConcreteInstance()){
    G4int n_trajectories = trajectoryContainer->entries();
    for(G4int i=0; i<n_trajectories; ++i){
      auto trj = (G4Trajectory*)((*(anEvent->GetTrajectoryContainer()))[i]);
      trj->DrawTrajectory();
    }
  }
#endif

  gAnaMan.EndOfEventAction();
}
