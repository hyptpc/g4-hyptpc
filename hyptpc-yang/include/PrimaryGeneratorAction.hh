/*
  PrimaryGeneratorAction.hh
  2017/10  Yang
*/

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include <vector>

class DetectorConstruction;
class AnalysisManager;
class EvtGen;
class EvtVector4R;
class G4ParticleGun;
class G4ParticleTable;
//class HepLorentzVector;
class EvtStdHep;
class EvtParticle;
class TFile;

typedef CLHEP::Hep3Vector ThreeVector;
typedef CLHEP::HepLorentzVector HepLorentzVector;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction( DetectorConstruction *det,
			  AnalysisManager *anaManager=0,
			  EvtGen *evtgen=0);
  ~PrimaryGeneratorAction();
  /*
  PrimaryGeneratorAction( DetectorConstruction *det,
			  Analysis *anaManager=0,
			  EvtGen *evtGen=0);
  ~PrimaryGeneratorAction();
  */
public:
  void GeneratePrimaries( G4Event *anEvent );
  void GenerateTest( G4Event *anEvent, G4ThreeVector D, G4ThreeVector P );
  void GenerateTest2( G4Event *anEvent, EvtGen *evt, G4ThreeVector D, G4ThreeVector P);
  void GenerateElasticPiM( G4Event *anEvent, G4ThreeVector D, G4ThreeVector P );
  void GenerateElasticPiP( G4Event *anEvent, G4ThreeVector D, G4ThreeVector P );
  void GenerateBeamPiM( G4Event *anEvent, G4ThreeVector D, G4ThreeVector P );
  void GenerateBeamPiP( G4Event *anEvent, G4ThreeVector D, G4ThreeVector P );
  void GeneratePiPSchannel( G4Event *anEvent, EvtGen *evt, G4ThreeVector D, G4ThreeVector P);
  void GenerateBeamKaonM( G4Event *anEvent, G4ThreeVector D, G4ThreeVector P );
  void GenerateBeamKaonMBeamRoot( G4Event *anEvent );
  void GenerateHDiBaryonBeamRoot_test( G4Event *anEvent );
  void GenerateHDiBaryon_test( G4Event *anEvent, EvtGen *evt, G4ThreeVector D, G4ThreeVector P);
  void GenerateBeamKaonMBr( G4Event *anEvent, G4ThreeVector D, G4ThreeVector P );
  void GenerateLambdaStarBr( G4Event *anEvent, EvtGen *evt, G4ThreeVector D, G4ThreeVector P);
  void GenerateLambdaStarBrCS( G4Event *anEvent, EvtGen *evt, G4ThreeVector D, G4ThreeVector P);

  double bw_swave( double beamp );
  void ReadBeamProfile( G4ThreeVector &X, G4ThreeVector &Y);

  G4ThreeVector GaussPosition( double sigx, double sigy, double hz );
  G4ThreeVector GaussMomentum( double sigx, double sigy, double sigz );
  G4ThreeVector GaussMomentum( double umean, double usig, double vmean, double vsig );
  double GaussZ( double hz );

  bool Scattering2Body( double Mi1, double Mi2, double Mf1, double Mf2, const G4ThreeVector & Pini, G4ThreeVector & Pf1, G4ThreeVector & Pf2, double cosx);

  bool Decay2Body( double Mini, double Mf1, double Mf2, const G4ThreeVector & Pini, G4ThreeVector & Pf1, G4ThreeVector & Pf2, double cosx );

  G4double Legendre(G4int order, G4double x);

  void GenerateScatteredKPlus(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P);

protected:
  DetectorConstruction *det_;
  AnalysisManager *anaMan_;
  EvtGen *evtgen_;
  //BeamParam *BP_;
  G4ParticleTable* particleTable;

  G4ParticleGun *chooseGun( int Pid );
  G4ParticleGun *chooseGun( const G4String & name );
  G4ParticleGun *particleGun;

  void makeGun(G4Event *anEvent, int partnum, EvtVector4R x, EvtVector4R y);
  void DeleteGuns();
  void GenerateDecay(G4Event* anEvent, EvtGen *evtGenerator, EvtParticle* particle, G4ThreeVector D);
  G4ThreeVector BeamMomRotate(G4ThreeVector beam, G4double x);

  // Beam //
  double bpx_; //beam momentum
  double bpy_;
  double bpz_;
  double bvx_; //beam vertex
  double bvy_;
  double bvz_;

  //TFile * beam_file;
  int bp_file_ndata;
  int bp_nAccess;


private:
  //AnalysisManager* anaMgr;
  
};

#endif
