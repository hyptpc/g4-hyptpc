// -*- C++ -*-

#ifndef TPCPADHELPER_HH
#define TPCPADHELPER_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

namespace TPCPadHelper
{
// Pad / Layer / Row indices are unified to start from 0
G4int GetPadID(const G4int layerID, const G4int rowID);
G4int GetLayerID(G4int padID);
G4int GetRowID(G4int padID);
G4double GetTheta(G4int padID);
G4double GetTheta(const G4int layerID, const G4double m_row);
G4double GetMrow(const G4int layerID, const G4double m_phi);
G4double GetRadius(const G4int layerID);
G4double GetLength(const G4int layerID);
G4double GetsTheta(const G4int layerID);
G4double GetR(G4int padID);
G4double GetZTarget();
G4ThreeVector GetPosition(G4int padID);
G4ThreeVector GetPosition(const G4int layerID, const G4double m_row);
// Find PadID from global position (z, x)
// Returns:
//    0 or positive : Valid PadID
//    -layer        : Hit inside the gap between layer and layer-1
//    -1000         : Not found (outside detector volume)
G4int FindPadID(G4double z, G4double x);
G4double ArcLength(const G4int layerID, const G4double row1, const G4double row2);
// Dead = pad on center frame / empty pad
G4bool IsDead(const G4int padID);
G4bool IsDead(const G4int layerID, const G4int rowID);
// Noise = abnormal waveform pad (E72)
G4bool Noise(const G4int padID);
G4bool Noise(const G4int layerID, const G4int rowID);
// Compatibility: always drop dead; noise only if drop_noise (TPCDropNoisePad).
G4bool GetDeadCon(const G4int padID, G4bool drop_noise = true);
G4bool GetDeadCon(const G4int layerID, const G4int rowID, G4bool drop_noise = true);

inline G4String ClassName()
{
  static G4String s_name("TPCPadHelper");
  return s_name;
}
}

#endif
