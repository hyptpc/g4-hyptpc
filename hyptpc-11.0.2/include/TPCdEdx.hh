// -*- C++ -*-

#ifndef TPCDEDX_HH
#define TPCDEDX_HH

#include <G4Types.hh>

namespace TPCdEdx
{
// Bethe-Bloch mean dE/dx for HypTPC P10 gas (ADC scale via TpcConversionFactor).
G4double dEdx(G4double mass /* G4: MeV/c2 */, G4double beta);

// Empirical dE/dx fluctuation width (pion/proton mass bands). Not Bethe-Bloch.
// Pass G4 internal units; converted to GeV inside for the empirical fit.
G4double dEdxSigma(G4double mass /* G4: MeV/c2 */, G4double mom /* G4: MeV/c */);
}

#endif
