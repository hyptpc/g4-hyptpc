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

  G4int legendre_order = 3;
  switch(pk_min_diff){
  case 723: //diff cross section = 0 below etaLambda threshold (exact value = 723.293)
    {
      G4double maximum_value = 0.;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 724:
    {
      G4double legendre_coeff[legendre_order] = {0.0162793, 0.00440176, 0.00677799};
      G4double maximum_value = 0.0274591;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 726:
    {
      G4double legendre_coeff[legendre_order] = {0.0297064, -0.000898291, 0.0138006};
      G4double maximum_value =  0.0444053;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 728:
    {
      G4double legendre_coeff[legendre_order] = {0.0445592, -0.00801745, 0.0139182};
      G4double maximum_value = 0.0664949;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 730:
    {
      G4double legendre_coeff[legendre_order] = {0.0657561, 0.0202863, 0.0526398};
      G4double maximum_value = 0.138682;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 732:
    {
      G4double legendre_coeff[legendre_order] = {0.0688638, -0.00504829, 0.0363991};
      G4double maximum_value = 0.110311;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 734:
    {
      G4double legendre_coeff[legendre_order] = {0.0902174, 0.0288109, 0.104857};
      G4double maximum_value = 0.223885;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
     }
  case 738:
    {
      G4double legendre_coeff[legendre_order] = {0.10436, 0.0235637, 0.0526217};
      G4double maximum_value = 0.180545;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 742:
    {
      G4double legendre_coeff[legendre_order] = {0.111125, -0.0107509, 0.0518198};
      G4double maximum_value = 0.173696;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 746:
    {
      G4double legendre_coeff[legendre_order] = {0.113412, -0.0159145, 0.0277801};
      G4double maximum_value = 0.157107;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 750:
    {
      G4double legendre_coeff[legendre_order] = {0.0944732, -0.000168073, 0.0425309};
      G4double maximum_value = 0.137172;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 754:
    {
      G4double legendre_coeff[legendre_order] = {0.101783, -0.0137408, 0.0162012};
      G4double maximum_value = 0.131725;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 758:
    {
      G4double legendre_coeff[legendre_order] = {0.0942948, -0.0109754, 0.0459009};
      G4double maximum_value = 0.151171;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 762:
    {
      G4double legendre_coeff[legendre_order] = {0.0723755, -0.044794, -0.00302102};
      G4double maximum_value = 0.114148;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 766:
    {
      G4double legendre_coeff[legendre_order] = {0.0787211, -0.0477353, 0.0156414};
      G4double maximum_value = 0.142098;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 770:
    {
      G4double legendre_coeff[legendre_order] = {0.0586327, 0.00726473, 0.0193116};
      G4double maximum_value = 0.0852091;  // maximum value of legendre func
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
