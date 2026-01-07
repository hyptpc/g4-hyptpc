// -*- C++ -*-

#ifndef KINEMATICS_HH
#define KINEMATICS_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include "TMinuit.h"

namespace Kinematics
{
  static const auto m12C = 11.174864*CLHEP::GeV;
  static const auto m11B = 10.252548*CLHEP::GeV;
  static const auto m10Be = 9.325504*CLHEP::GeV;
  G4double DecayMomentum(double M,double m1, double m2);
  G4ThreeVector HarmonicFermiMomentum( G4int type );
  G4int         HarmonicFermiMomentumDeuteron( G4double* Kf );
  G4double      Legendre( G4int order, G4double x );
  G4ThreeVector SphericalRandom( void );
  G4ThreeVector FermiGasMomentum(double p_f = 250);//MeV
  G4ThreeVector FermiEmphiricalMomentum();//MeV
  void FermiScattering(TLorentzVector p_in, TLorentzVector& p_out, TLorentzVector& p_proton, 
      int offshellness = 0,int conf = 0,
      double mass_nuclei = m11B, double mass_daughter = m10Be);//MeV
	G4ThreeVector MultitrackVertex(int ntrack, double *x0, double *y0,
			    double *u0, double *v0);

	G4ThreeVector RotateAlongBeam(G4ThreeVector Beam, G4ThreeVector Vect, G4double phi);
  inline G4String ClassName( void )
  {
    static G4String s_name("Kinematics");
    return s_name;
  }
}

#endif
