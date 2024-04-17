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
G4bool SigmaZeroPiZeroBC(const G4double cos_theta, const G4double pk);
G4bool SigmaZeroPiZeroCB(const G4double cos_theta, const G4double pk);
G4bool SigmaMinusPiPlus(const G4double cos_theta, const G4double pk);
G4bool SigmaPlusPiMinus(const G4double cos_theta, const G4double pk);
G4bool KpElastic(const G4double cos_theta, const G4double pk);

inline G4String ClassName()
{
  static G4String s_name("DiffCrossSection");
  return s_name;
}

}

#endif
