// -*- C++ -*-

#ifndef DIFFCROSSSECTION_HH
#define DIFFCROSSSECTION_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

namespace DiffCrossSection
{
G4bool LambdaEta(const G4double cos_theta, const G4double pk);
G4bool LambdaPiZeroBC(const G4double cos_theta, const G4double pk);
G4bool LambdaPiZeroCB(const G4double cos_theta, const G4double pk);
G4bool SigmaZeroPiZero(const G4double cos_theta, const G4double pk);

inline G4String ClassName()
{
  static G4String s_name("DiffCrossSection");
  return s_name;
}

}

#endif
