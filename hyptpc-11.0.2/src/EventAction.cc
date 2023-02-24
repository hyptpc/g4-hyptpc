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
#include "BACSD.hh"
#include "ConfMan.hh"
#include "FuncName.hh"
#include "BH2SD.hh"
#include "FTOFSD.hh"
#include "HTOFSD.hh"
#include "LACSD.hh"
#include "TPCSD.hh"
#include "SCHSD.hh"
#include "SDCSD.hh"
#include "TargetSD.hh"
#include "VPSD.hh"
#include "WCSD.hh"

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
  // G4cout << FUNC_NAME << G4endl;
  gAnaMan.BeginOfEventAction();
}

//_____________________________________________________________________________
void
EventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4int eventID = anEvent-> GetEventID();
  if(eventID % 100 == 0){
    G4cout << FUNC_NAME << G4endl
	   << "   Event number = " << eventID << G4endl;
  }
  G4SDManager* SDManager= G4SDManager::GetSDMpointer();

  // get "Hit Collection of This Event"
  G4HCofThisEvent* HCTE= anEvent-> GetHCofThisEvent();
  if(! HCTE) return;  // no hits in this events. nothing to do!

  static const G4int id_bh2 = SDManager-> GetCollectionID("BH2/hit");
  if(id_bh2 >= 0){
    const auto HC = (G4THitsCollection<BH2Hit>*)HCTE->GetHC(id_bh2);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetBH2Data((*HC)[i]);
    }
  }

  static const G4int idcounter = SDManager->GetCollectionID("TPC/hit");
  if(idcounter > 0){
    auto hc = (G4THitsCollection<TPCHit>*)HCTE->GetHC(idcounter);
    G4int nhits= hc -> entries();
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
    for(G4int i=0; i<nhits; ++i){
      G4ThreeVector vtxpos = (*hc)[i]-> GetVtxPosition();
      G4ThreeVector vtxmom = (*hc)[i]-> GetVtxMomentum();
      G4double vtxene =(*hc)[i]-> GetVtxEnergy();
      G4ThreeVector xyz = (*hc)[i]-> GetPosition();
      G4ThreeVector mom = (*hc)[i]-> GetMomentum();
      G4double tof= (*hc)[i]-> GetTOF();
      G4int tid = (*hc)[i]-> GetTrackID();
      G4int ptid = (*hc)[i]-> GetParentID();
      G4int ptid_pid = (*hc)[i]-> GetParentID_pid();
      G4int pid = (*hc)[i]-> GetParticleID();
      G4double mass = (*hc)[i]-> GetMass();
      G4int charge = (*hc)[i]-> GetCharge();
      // std::cout<<"pid="<<pid<<", mass="<<mass<<", charge="<<charge<<std::endl;
      // getchar();
      G4int ilay = (*hc)[i]-> GetPadLay();
      //      G4double mass = (*hc)[i]-> GetPDGMass(); //mass(GeV)
      G4int parentid = (*hc)[i]-> GetParentID();
      G4double tlength = (*hc)[i]-> GettLength();
      G4int irow=(*hc)[i]-> GetPadRow();
      G4double beta = (*hc)[i]-> GetBeta();
      G4double edep = (*hc)[i]-> GetEdep();
      G4double slength = (*hc)[i]-> GetsLength();
      //      G4VTrajectoryPoint *tp -> hc->GetPoint(i);
      //    G4int nhits= hc -> entries();
      if(nparticle==0){
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

      }else if(nparticle>0){
	//	G4cout<<nparticle<<G4endl;
      }
      if((pidtr[nparticle-1] != pid) || (pidtr[nparticle-1] == pid && vtxpxtpc[nparticle-1] != vtxmom[0] && vtxpytpc[nparticle-1] != vtxmom[1] && vtxpztpc[nparticle-1] != vtxmom[2])){
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
        if(ptidtpc[nparticle-1] != ptid){
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
	gAnaMan.SetCounterData(nparticle-1,tof, xyz, mom, tid, pid, ilay,
                               irow, beta, edep/CLHEP::MeV, parentid,
                               tlength,slength);
      }
    }
    for(G4int i=0;i<nparticle;i++){
      gAnaMan.SetTPCData(i, pidtr[i], ptidtpc[i], ptidtpc_pid[i], pxtpc[i],pytpc[i],pztpc[i],pptpc[i], qqtpc[i], pmtpc[i], detpc[i], lentpc[i],laytpc[i],
                         vtxpxtpc[i],vtxpytpc[i],vtxpztpc[i],
                         vtxxtpc[i],vtxytpc[i],vtxztpc[i], vtxenetpc[i]);
    }
  }

  static const G4int id_target = SDManager->GetCollectionID("TGT/hit");
  if(id_target >= 0){
    const auto HC = (G4THitsCollection<TargetHit>*)HCTE->GetHC(id_target);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetTargetData((*HC)[i]);
    }
  }

  static const G4int id_htof = SDManager-> GetCollectionID("HTOF/hit");
  if(id_htof >= 0){
    const auto HC = (G4THitsCollection<HTOFHit>*)HCTE->GetHC(id_htof);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetHTOFData((*HC)[i]);
    }
  }

  // if(gConf.Get<G4int>("UseAC") == 1.){
  //   static const G4int id_ac = SDManager->GetCollectionID("AC/hit");
  //   if(id_ac >= 0){
  //     const auto HC = (G4THitsCollection<BACHit>*)HCTE->GetHC(id_ac);
  //     for(G4int i=0, n=HC->entries(); i<n; ++i){
  // 	gAnaMan.SetACData((*HC)[i]);
  //     }
  //   }
  // }

  static const G4int id_sdc = SDManager-> GetCollectionID("SDC/hit");
  if(id_sdc >= 0){
    const auto HC = (G4THitsCollection<SDCHit>*)HCTE->GetHC(id_sdc);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetSDCData((*HC)[i]);
    }
  }

  static const G4int id_sch = SDManager-> GetCollectionID("SCH/hit");
  if(id_sch >= 0){
    const auto HC = (G4THitsCollection<SCHHit>*)HCTE->GetHC(id_sch);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetSCHData((*HC)[i]);
    }
  }

  static const G4int id_ftof = SDManager-> GetCollectionID("FTOF/hit");
  if(id_ftof >= 0){
    const auto HC = (G4THitsCollection<FTOFHit>*)HCTE->GetHC(id_ftof);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetFTOFData((*HC)[i]);
    }
  }

  static const G4int id_lac = SDManager-> GetCollectionID("LAC/hit");
  if(id_lac >= 0){
    const auto HC = (G4THitsCollection<LACHit>*)HCTE->GetHC(id_lac);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetLACData((*HC)[i]);
    }
  }

  static const G4int id_wc = SDManager-> GetCollectionID("WC/hit");
  if(id_wc >= 0){
    const auto HC = (G4THitsCollection<WCHit>*)HCTE->GetHC(id_wc);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetWCData((*HC)[i]);
    }
  }

  static const G4int id_vp = SDManager-> GetCollectionID("VP/hit");
  if(id_vp >= 0){
    const auto HC = (G4THitsCollection<VPHit>*)HCTE->GetHC(id_vp);
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      gAnaMan.SetVPData((*HC)[i]);
    }
  }

  auto trajectoryContainer = anEvent->GetTrajectoryContainer();
  if(trajectoryContainer && G4VVisManager::GetConcreteInstance()){
    G4int n_trajectories = trajectoryContainer->entries();
    for(G4int i=0; i<n_trajectories; ++i){
      auto trj = (G4Trajectory*)((*(anEvent->GetTrajectoryContainer()))[i]);
      trj->DrawTrajectory();
    }
  }
  gAnaMan.EndOfEventAction();
}
