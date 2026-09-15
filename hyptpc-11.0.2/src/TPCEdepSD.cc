// -*- C++ -*-

#include "TPCEdepSD.hh"

#include <G4VPhysicalVolume.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VTouchable.hh>
#include <G4TouchableHistory.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4DynamicParticle.hh>
#include <G4DecayProducts.hh>
#include <G4PhysicsLogVector.hh>
#include <G4ParticleChangeForDecay.hh>
#include <G4DecayProcessType.hh>
#include <Randomize.hh>

#include <cmath>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "TPCEdepHit.hh"
#include "TPCPadHelper.hh"



namespace
{
  const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
TPCEdepSD::TPCEdepSD( const G4String& name )
  : G4VSensitiveDetector( name ),
    m_hits_collection()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
TPCEdepSD::~TPCEdepSD( void )
{
}

//_____________________________________________________________________________
void
TPCEdepSD::Initialize( G4HCofThisEvent* HCTE )
{
  ntrk=0;
  m_hits_collection = new G4THitsCollection<TPCEdepHit>( SensitiveDetectorName,
						     collectionName[0] );
  // push H.C. to "Hit Collection of This Event"
  G4int hcid = GetCollectionID(0);
  HCTE->AddHitsCollection( hcid, m_hits_collection );
}

//_____________________________________________________________________________
G4bool
TPCEdepSD::ProcessHits( G4Step* aStep, G4TouchableHistory* /* ROhist */ )
{

  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  const G4StepPoint* postStepPoint= aStep-> GetPostStepPoint();
  const G4Track* aTrack = aStep->GetTrack();
  G4int copyNo = preStepPoint -> GetPhysicalVolume()->GetCopyNo();
  //std::cout<<"copyNo : "<<copyNo<<std::endl;
  G4int copyNo_post = postStepPoint ->GetPhysicalVolume()->GetCopyNo();
  //if(copyNo!=copyNo_post) return false;

  
  //if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;
  //  if(preStepPoint-> GetStepStatus() == fGeomBoundary){
  G4String particleName;
  
	if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
    return false;
	
  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();

  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4ThreeVector pos= preStepPoint-> GetPosition();
  G4double hitx=pos.getX();
  G4double hitz=pos.getZ();

    G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack -> GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  G4ThreeVector mom= preStepPoint-> GetMomentum();
  G4double tof= preStepPoint-> GetGlobalTime();
  G4double beta= preStepPoint-> GetBeta();
  G4int tid =  aStep-> GetTrack()-> GetTrackID();
  G4int pid =  aStep-> GetTrack()-> GetDefinition() -> GetPDGEncoding();
  //  if(abs(pid)==11)return false; //remove e+, e-
  
  G4double mass = aStep -> GetTrack()->GetDynamicParticle()->GetMass();
  //  G4cout<<mass<<G4endl;
  G4int charge = aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge();
  //  G4double edep = aStep->GetTotalEnergyDeposit();
  //  G4double edep = aStep->GetNonIonizingEnergyDeposit();
  //  G4double deltaep = aStep->GetDeltaEnergy();
  G4double tlength = aStep->GetTrack()-> GetTrackLength();
  G4double slength = aStep->GetTrack()-> GetStepLength();
  // G4double length = aStep-> GetStepLength();
  G4int parentID =  aStep-> GetTrack()-> GetParentID();

  G4int parentID_pid = aStep-> GetTrack()->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();

  G4int iLay_copyNo=copyNo - 2000;
  //G4int iPad = padHelper::findPadID(hitz, hitx);
  //G4int iLay= padHelper::getLayerID(iPad);
  G4int iPad = TPCPadHelper::FindPadID(hitz, hitx);
  // valid layers 0..31; GetLayerID sentinel 32 = neg / OOB padID
  if (TPCPadHelper::GetLayerID(iPad) > 31)
    return false;
  //G4int iLay = TPCPadHelper::GetLayerID(iPad);
  G4int iLay = iLay_copyNo;
  //G4int iRow= padHelper::getRowID(iPad);
  G4int iRow = TPCPadHelper::GetRowID(iPad);
  
	G4double PadLen = TPCPadHelper::GetLength(iLay);
	G4ThreeVector PadPos(hitx,0,hitz - TPCPadHelper::GetZTarget()); 
	G4ThreeVector MomT(mom.x(),0,mom.z());	
	G4double alpha = PadPos.theta()-MomT.theta();
	G4double PathT = PadLen * 1./cos(alpha);
	G4double Pitch = mom.y()/MomT.mag();
	G4double Path = PathT * sqrt(1+Pitch*Pitch);
	slength = Path;

	/*
	G4double edepMean =TPCdEdx(mass,beta)*Path; 
	//	G4double IonEn = (0.9 * 188 + 0.1 * 41.7)*eV;
//	G4double edepSig = sqrt(edepMean / IonEn ) * IonEn;
	G4double edepSig =TPCdEdxSig(mass,mom.mag())*Path;
	if(edepSig / edepMean < 0.01 or edepSig/edepMean > 0.5)edepSig = 0.2*edepMean;
	G4double edep = G4RandGauss::shoot(edepMean,edepSig);
	if(edep <0.1* edepMean)edep = 0.1*edepMean;
	*/

	G4double conversion_factor = 7388.11; //MeV/cm -> ADC/mm 
	G4double cmTomm = 10;
	G4double edep = aStep->GetTotalEnergyDeposit() * conversion_factor * cmTomm; 


	//if(iLay_copyNo==4 && pos.getY()>-1 * 97./2.)return false;
	
#ifdef DEBUG

  // Compare geometry layer vs FindPadID-derived layer (not self-compare of iLay).
  const G4int iLay_pad = TPCPadHelper::GetLayerID(iPad);
  const G4double radius = std::hypot(hitx, hitz - TPCPadHelper::GetZTarget());

  G4ThreeVector pos_post = postStepPoint->GetPosition();
  const G4double hitx_post = pos_post.getX();
  const G4double hitz_post = pos_post.getZ();
  const G4int iPad_post = TPCPadHelper::FindPadID(hitz_post, hitx_post);
  const G4int iLay_post = TPCPadHelper::GetLayerID(iPad_post);
  const G4int iLay_post_copy = copyNo_post - 2000;

  if (iLay_pad < 0 || iLay_pad > 31 || iLay_pad != iLay_copyNo) {
    G4cout << "pre strange layer"
           << " radius=" << radius
           << " iPad=" << iPad
           << " iLay_copyNo=" << iLay_copyNo
           << " iLay_pad=" << iLay_pad
           << " hitx=" << hitx
           << " hitz=" << hitz
           << G4endl;
  }
  // Post: only when still in TPCEdep copy range, or FindPadID returned a sentinel.
  const G4bool post_in_edep =
    (copyNo_post >= 2000 && copyNo_post < 2000 + 32);
  if (iPad_post < 0
      || (post_in_edep
          && (iLay_post < 0 || iLay_post > 31 || iLay_post != iLay_post_copy))) {
    G4cout << "post strange layer"
           << " radius=" << std::hypot(hitx_post, hitz_post - TPCPadHelper::GetZTarget())
           << " iPad=" << iPad_post
           << " iLay_copyNo=" << iLay_post_copy
           << " iLay_pad=" << iLay_post
           << " hitx=" << hitx_post
           << " hitz=" << hitz_post
           << G4endl;
  }

#endif


  
  G4String name = physVol->GetName();
  
  



  if(name=="TPC_PV"){
    name="EdepPV-1";
  }

  sscanf(name,"EdepPV%d",&iLay_copyNo);
  TPCEdepHit* ahit= new TPCEdepHit(pos, mom, tof, tid, pid, iLay, iRow, beta, edep,parentID,tlength, mass, charge, VertexPosition, VertexMomentum, VertexEnergy,slength, parentID_pid );

  m_hits_collection-> insert(ahit);
  return true;

}

//_____________________________________________________________________________
void
TPCEdepSD::EndOfEvent(G4HCofThisEvent* /* HCTE */ )
{
}

//_____________________________________________________________________________
void
TPCEdepSD::DrawAll( void )
{
}

G4double
TPCEdepSD::TPCdEdx(G4double mass/*MeV/c2*/, G4double beta){

  G4double rho=0.; //[g cm-3]
  G4double ZoverA=0.; //[mol g-1]
  G4double I=0.; //[eV]
  G4double density_effect_par[6]={0.}; //Sternheimer’s parameterization
  //P10  
	rho = std::pow(10.,-3)*(0.9*1.662 + 0.1*0.6672);
	ZoverA = 17.2/37.6;
	I = 0.9*188.0 + 0.1*41.7;
	density_effect_par[0] = 0.9*0.19714 + 0.1*0.09253;
	density_effect_par[1] = 0.9*2.9618 + 0.1*3.6257;
	density_effect_par[2] = 0.9*1.7635 + 0.1*1.6263;
	density_effect_par[3] = 0.9*4.4855 + 0.1*3.9716;
	density_effect_par[4] = 0.9*11.9480 + 0.1*9.5243;
	density_effect_par[5] = 0.;

  G4double Z = 1.;
  G4double me = 0.5109989461; //[MeV]
  G4double K = 0.307075; //[MeV cm2 mol-1]
  G4double constant = rho*K*ZoverA; //[MeV cm-1]
  constant = constant;
  G4double I2 = I*I; //Mean excitaion energy [eV]
  G4double beta2 = beta*beta;
  G4double gamma2 = 1./(1.-beta2);
  G4double MeVToeV = std::pow(10.,6);
  G4double Wmax = 2*me*beta2*gamma2/((me/mass+1.)*(me/mass+1.)+2*(me/mass)*(std::sqrt(gamma2)-1));
  G4double delta = DensityEffectCorrection(std::sqrt(beta2*gamma2), density_effect_par);
  G4double dedx = constant*Z*Z/beta2*(0.5*std::log(2*me*beta2*gamma2*Wmax*MeVToeV*MeVToeV/I2) - beta2 - 0.5*delta);
  
	
  G4double conversion_factor = 11073.3;
  return conversion_factor*dedx;

}
G4double
TPCEdepSD::TPCdEdxSig(G4double mass/*MeV/c2*/, G4double mom){
	G4double par[3]={0,0,0} ;
	if(mass < 0.2){//pion;
		par[0] = 7.792;
		par[1] =	-8.704;
		par[2] = 4.477;
	}
	else if(mass > 0.7){//proton;
		par[0] = 33.92;
		par[1] = -26.24;
		par[2] = 6.259; 
	}
	G4double value = par[0]+par[1]*mom+par[2]*mom*mom;
	return value;
}
G4double
TPCEdepSD::DensityEffectCorrection(G4double betagamma, G4double *par){

    //reference : Sternheimer’s parameterizatio(PDG)
    //notation : par[0] : a, par[1] : k, par[2] : x0, par[3] : x1, par[4] : _C, par[5] : delta0
    G4double constant = 2*std::log(10.);
    G4double delta = 0.;
    G4double X = std::log10(betagamma);
    if(X<=par[2]) delta = par[5]*std::pow(10., 2*(X - par[2]));
    else if(par[2]<X && X<par[3]) delta = constant*X - par[4] + par[0]*std::pow((par[3] - X), par[1]);
    else if(X>=par[3]) delta = constant*X - par[4];

  return delta;

}
//_____________________________________________________________________________
void
TPCEdepSD::PrintAll( void )
{
  m_hits_collection-> PrintAllHits();
}

