// -*- C++ -*-

#ifndef KINEMATICS_HH
#define KINEMATICS_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

namespace Kinematics
{
G4ThreeVector HarmonicFermiMomentum(G4int type);
G4int         HarmonicFermiMomentumDeuteron(G4double* Kf);
G4double      Legendre(G4int order, G4double x);
G4double      EffectiveThickness(const G4ThreeVector pos, const G4ThreeVector mom, const G4double target_center_z, const G4double target_r);
G4ThreeVector RandomVertex(const G4ThreeVector pos, const G4ThreeVector mom, const G4double target_center_z, const G4double target_r, const G4double target_height);

inline G4String ClassName()
{
  static G4String s_name("Kinematics");
  return s_name;
}
}

#endif
