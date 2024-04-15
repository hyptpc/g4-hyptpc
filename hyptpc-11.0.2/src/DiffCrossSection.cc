// -*- C++ -*-

#include "DiffCrossSection.hh"

#include <Randomize.hh>

#include "AnaManager.hh"
#include "FuncName.hh"
#include "Kinematics.hh"

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
}

namespace DiffCrossSection
{

//______________________________________________________________________________
// calculate diff. cross sec. using CB data (DOI: 10.1103/PhysRevC.64.055205) 
// cross section  = 0 before eta, lambda threshold
G4bool
LambdaEta(const G4double cos_theta, const G4double pk)
{
  std::vector<G4int> pk_list{723, 724, 726, 728, 730, 732, 734, 738, 742, 746, 750, 754, 758, 762, 766, 770};
  G4int pk_min_diff = pk_list[0];
  G4double min_diff = 1000;  // initialize
  for (G4int i = 0, n = pk_list.size(); i < n; i++) {
    G4double diff = std::abs( pk-pk_list[i] );
    if (diff < min_diff) {
      min_diff    = diff;
      pk_min_diff = pk_list[i];
    }
  }
  gAnaMan.SetDiffCrossSecMom(pk_min_diff);

  G4int legendre_order = 3;  // 3 means using up to 2nd order (there are 3 terms)
  switch(pk_min_diff){
  case 723: //diff cross section = 0 below etaLambda threshold (exact value = 723.293)
    {
      G4double maximum_value = 0.;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 724:
    {
      G4double legendre_coeff[legendre_order] = {0.0162205, 0.00440677, 0.00685386};
      G4double maximum_value = 0.0274811;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 726:
    {
      G4double legendre_coeff[legendre_order] = {0.0295538, -0.000886306, 0.0139095};
      G4double maximum_value = 0.0443497;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 728:
    {
      G4double legendre_coeff[legendre_order] = {0.0444324, -0.00800496, 0.0140956};
      G4double maximum_value = 0.066533;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 730:
    {
      G4double legendre_coeff[legendre_order] = {0.0652961, 0.0202536, 0.0534762};
      G4double maximum_value = 0.139026;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 732:
    {
      G4double legendre_coeff[legendre_order] = {0.0685303, -0.00503868, 0.0368818};
      G4double maximum_value = 0.110451;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 734:
    {
      G4double legendre_coeff[legendre_order] = {0.0893263, 0.0287288, 0.106272};
      G4double maximum_value = 0.224327;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
     }
  case 738:
    {
      G4double legendre_coeff[legendre_order] = {0.103881, 0.0235218, 0.0532157};
      G4double maximum_value = 0.180619;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 742:
    {
      G4double legendre_coeff[legendre_order] = {0.110574, -0.0107958, 0.0524904};
      G4double maximum_value = 0.173861;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 746:
    {
      G4double legendre_coeff[legendre_order] = {0.113162, -0.0159179, 0.0281021};
      G4double maximum_value = 0.157182;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 750:
    {
      G4double legendre_coeff[legendre_order] = {0.094048, -0.000138357, 0.0430225};
      G4double maximum_value = 0.137209;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 754:
    {
      G4double legendre_coeff[legendre_order] = {0.101632, -0.0137436, 0.0165167};
      G4double maximum_value = 0.131892;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 758:
    {
      G4double legendre_coeff[legendre_order] = {0.0938071, -0.0109858, 0.0463662};
      G4double maximum_value = 0.151159;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 762:
    {
      G4double legendre_coeff[legendre_order] = {0.0724061, -0.0447803, -0.00309023};
      G4double maximum_value = 0.114096;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 766:
    {
      G4double legendre_coeff[legendre_order] = {0.0785833, -0.0477332, 0.0158581};
      G4double maximum_value = 0.142175;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 770:
    {
      G4double legendre_coeff[legendre_order] = {0.0584505, 0.00726072, 0.0196413};
      G4double maximum_value = 0.0853524;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }
}

}
