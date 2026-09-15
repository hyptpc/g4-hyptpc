// -*- C++ -*-

#include "TPCdEdx.hh"

#include <cmath>

#include "ConfMan.hh"
#include "TPCPadHelper.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();

// Sternheimer density-effect correction (PDG parameterization).
// par[0]=a, [1]=k, [2]=x0, [3]=x1, [4]=C, [5]=delta0
G4double
DensityEffectCorrection(G4double betagamma, G4double* par)
{
  const G4double constant = 2.0 * std::log(10.0);
  G4double delta = 0.0;
  const G4double X = std::log10(betagamma);
  if (X <= par[2]) {
    delta = par[5] * std::pow(10.0, 2.0 * (X - par[2]));
  }
  else if (par[2] < X && X < par[3]) {
    delta = constant * X - par[4]
            + par[0] * std::pow((par[3] - X), par[1]);
  }
  else if (X >= par[3]) {
    delta = constant * X - par[4];
  }
  return delta;
}
}

namespace TPCdEdx
{

//_____________________________________________________________________________
// Mean dE/dx for HypTPC P10 gas (Ar 90% + CH4 10%), ADC scale via
// TpcConversionFactor. Same formula as analyzer Kinematics::HypTPCdEdx(0,...).
G4double
dEdx(G4double mass /*MeV/c2*/, G4double beta)
{
  // --- P10 gas parameters (weighted Ar:CH4 = 0.9:0.1) ---
  const G4double rho =
    std::pow(10.0, -3.0) * (0.9 * 1.662 + 0.1 * 0.6672); // [g cm-3]
  const G4double ZoverA = 17.2 / 37.6; // [mol g-1]
  const G4double I = 0.9 * 188.0 + 0.1 * 41.7; // [eV]
  G4double density_effect_par[6] = {0.0};
  density_effect_par[0] = 0.9 * 0.19714 + 0.1 * 0.09253;
  density_effect_par[1] = 0.9 * 2.9618 + 0.1 * 3.6257;
  density_effect_par[2] = 0.9 * 1.7635 + 0.1 * 1.6263;
  density_effect_par[3] = 0.9 * 4.4855 + 0.1 * 3.9716;
  density_effect_par[4] = 0.9 * 11.9480 + 0.1 * 9.5243;
  density_effect_par[5] = 0.0;

  const G4double Z = 1.0;
  const G4double me = 0.5109989461; // [MeV]
  const G4double K = 0.307075; // [MeV cm2 mol-1]
  const G4double constant = rho * K * ZoverA; // [MeV cm-1]
  const G4double I2 = I * I;
  const G4double beta2 = beta * beta;
  const G4double gamma2 = 1.0 / (1.0 - beta2);
  const G4double MeVToeV = std::pow(10.0, 6.0);

  // Maximum energy transfer to an electron in a single collision
  const G4double me_over_m = me / mass;
  const G4double wmax_num = 2.0 * me * beta2 * gamma2;
  const G4double wmax_den =
    (me_over_m + 1.0) * (me_over_m + 1.0)
    + 2.0 * me_over_m * (std::sqrt(gamma2) - 1.0);
  const G4double Wmax = wmax_num / wmax_den;

  const G4double betagamma = std::sqrt(beta2 * gamma2);
  const G4double delta = DensityEffectCorrection(betagamma, density_effect_par);

  const G4double log_arg = 2.0 * me * beta2 * gamma2 * Wmax * MeVToeV * MeVToeV / I2;
  const G4double dedx = constant * Z * Z / beta2 * (0.5 * std::log(log_arg) - beta2 - 0.5 * delta);

  static const G4double conversion_factor =
    gConf.GetOrDefault<G4double>("TpcConversionFactor",
                                 TPCPadHelper::kDefaultConversionFactor);
  return conversion_factor * dedx;
}

//_____________________________________________________________________________
G4double
dEdxSigma(G4double mass /*MeV/c2*/, G4double mom)
{
  G4double par[3] = {0.0, 0.0, 0.0};
  if (mass < 0.2) { // pion
    par[0] = 7.792;
    par[1] = -8.704;
    par[2] = 4.477;
  }
  else if (mass > 0.7) { // proton
    par[0] = 33.92;
    par[1] = -26.24;
    par[2] = 6.259;
  }
  return par[0] + par[1] * mom + par[2] * mom * mom;
}

}
