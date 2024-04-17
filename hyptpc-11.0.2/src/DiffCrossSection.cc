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
// coefficients are determined by fitting ourselvs
// cross section = 0 before eta, lambda threshold
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

  G4double maximum_value = 0.;  // maximum value of legendre func
  std::vector<G4double> legendre_coeff;
  legendre_coeff.clear();
  
  switch(pk_min_diff){
  case 723: //diff cross section = 0 below etaLambda threshold (exact value = 723.293)
    maximum_value = 0.;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.});
    break;
  case 724:
    maximum_value = 0.0274811;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0162205, 0.00440677, 0.00685386});
    break;
  case 726:
    maximum_value = 0.0443497;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0295538, -0.000886306, 0.0139095});
    break;
  case 728:
    maximum_value = 0.066533;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0444324, -0.00800496, 0.0140956});
    break;
  case 730:
    maximum_value = 0.139026;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0652961, 0.0202536, 0.0534762});
    break;
  case 732:
    maximum_value = 0.110451;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0685303, -0.00503868, 0.0368818});
    break;
  case 734:
    maximum_value = 0.224327;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0893263, 0.0287288, 0.106272});
    break;
  case 738:
    maximum_value = 0.180619;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.103881, 0.0235218, 0.0532157});
    break;
  case 742:
    maximum_value = 0.173861;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.110574, -0.0107958, 0.0524904});
    break;
  case 746:
    maximum_value = 0.157182;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.113162, -0.0159179, 0.0281021});
    break;
  case 750:
    maximum_value = 0.137209;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.094048, -0.000138357, 0.0430225});
    break;
  case 754:
    maximum_value = 0.131892;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.101632, -0.0137436, 0.0165167});
    break;
  case 758:
    maximum_value = 0.151159;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0938071, -0.0109858, 0.0463662});
    break;
  case 762:
    maximum_value = 0.114096;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0724061, -0.0447803, -0.00309023});
    break;
  case 766:
    maximum_value = 0.142175;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0785833, -0.0477332, 0.0158581});
    break;
  case 770:
    maximum_value = 0.0853524;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.0584505, 0.00726072, 0.0196413});
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }

  G4double legendre_cos_theta = 0.;
  for (G4int order = 0, n_order = legendre_coeff.size(); order < n_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
  return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;

} //LambdaEta

//______________________________________________________________________________
// calculate diff. cross sec. using of Bubble Chamber data (https://doi.org/10.1016/0550-3213(70)90461-X)
// coefficients are determined by fitting ourselvs
G4bool
LambdaPiZeroBC(const G4double cos_theta, const G4double pk)
{
  std::vector<G4int> pk_list{617, 637, 658, 677, 699, 719, 740, 761, 773, 793};
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

  G4double maximum_value = 0.;  // maximum value of legendre func
  std::vector<G4double> legendre_coeff;
  legendre_coeff.clear();
  switch(pk_min_diff){
  case 617:
    maximum_value = 0.385164;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.184747, 0.123287, 0.136187, -0.0287959, -0.0302612});
    break;
  case 637:
    maximum_value = 0.587236;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.211184, 0.132109, 0.198364, 0.10548, -0.0599012});
    break;
  case 658:
    maximum_value = 0.583402;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.201419, 0.0922841, 0.198598, 0.0938456, -0.00274454});
    break;
  case 677:
    maximum_value = 0.608879;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.211806, 0.0539344, 0.195027, 0.171799, -0.0236871});
    break;
  case 699:
    maximum_value = 0.517189;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.176029, 0.0870065, 0.203049, 0.0983589, -0.0472533});
    break;
  case 719:
    maximum_value = 0.876051;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.236829, 0.11063, 0.378043, 0.132433, 0.0181166});
    break;
  case 740:
    maximum_value = 0.953612;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.229397, 0.12905, 0.342048, 0.151436, 0.101681});
    break;
  case 761:
    maximum_value = 0.821861;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.211328, 0.0863754, 0.422316, 0.0292175, 0.0726246});
    break;
  case 773:
    maximum_value = 1.00163;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.236279, 0.0850016, 0.461273, 0.113247, 0.105835});
    break;
  case 793:
    maximum_value = 1.15354;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.238692, 0.119329, 0.511488, 0.164256, 0.119778});
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }
  G4double legendre_cos_theta = 0.;
  for (G4int order = 0, n_order = legendre_coeff.size(); order < n_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
  return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;

} //LambdaPiZeroBC

//______________________________________________________________________________
// calculate diff. cross sec. using of Crystal Ball data (DOI: 10.1103/PhysRevC.80.025204)
// coefficients are used their results
G4bool
LambdaPiZeroCB(const G4double cos_theta, const G4double pk)
{
  std::vector<G4int> pk_list{514, 560, 581, 629, 659, 687, 714, 750};
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

  G4double maximum_value = 0.;  // maximum value of legendre func
  std::vector<G4double> legendre_coeff;
  legendre_coeff.clear();
  switch(pk_min_diff){
  case 514:
    maximum_value = 0.749;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.2240, 0.2740, 0.1800, 0.0750, -0.0040});
    break;
  case 560:
    maximum_value = 0.6005;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.2195, 0.2080, 0.1480, 0.0560, -0.0310});
    break;
  case 581:
    maximum_value = 0.511;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.2160, 0.1780, 0.1200, 0.0440, -0.0470});
    break;
  case 629:
    maximum_value = 0.5718;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.2171, 0.1367, 0.1720, 0.0980, -0.0520});
    break;
  case 659:
    maximum_value = 0.7811;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.2291, 0.1380, 0.2080, 0.1610, 0.0040, 0.0410});
    break;
  case 687:
    maximum_value = 0.7631;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.2411, 0.0980, 0.2350, 0.1720, -0.0040, 0.0210});
    break;
  case 714:
    maximum_value = 0.9277;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.2527, 0.0940, 0.3100, 0.1880, 0.0410, 0.0420});
    break;
  case 750:
    maximum_value = 0.9066;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.2619, 0.0670, 0.3830, 0.1350, 0.0600, -0.0003});
    break;
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }
  G4double legendre_cos_theta = 0.;
  for (G4int order = 0, n_order = legendre_coeff.size(); order < n_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
  return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;

} //LambdaPiZeroCB

//______________________________________________________________________________
// calculate diff. cross sec. using of Bubble Chamber data (https://doi.org/10.1016/0550-3213(70)90461-X)
// coefficients are determined by fitting ourselvs
G4bool
SigmaZeroPiZeroBC(const G4double cos_theta, const G4double pk)
{
  std::vector<G4int> pk_list{617, 637, 658, 677, 699, 719, 740, 761, 773, 793};
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

  G4double maximum_value = 0.;  // maximum value of legendre func
  std::vector<G4double> legendre_coeff;
  legendre_coeff.clear();
  switch(pk_min_diff){
  case 617:
    maximum_value = 0.357797;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.157887, -0.087196, 0.0797751, 0.0332213, 0.06616});
    break;
  case 637:
    maximum_value = 0.263738;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.111551, -0.0426781, 0.0689817, -0.0168688, 0.0236574});
    break;
  case 658:
    maximum_value = 0.380357;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.138992, -0.0583661, 0.125422, -0.0534204, 0.00415687});
    break;
  case 677:
    maximum_value = 0.372365;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.14989, -0.0701137, 0.0929754, -0.102943, -0.0435574});
    break;
  case 699:
    maximum_value = 0.3718;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.175131, -0.0356463, 0.133603, -0.0673572, -0.0399371});
    break;
  case 719:
    maximum_value = 0.549651;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.178271, -0.0203391, 0.245218, -0.0383859, 0.0674362});
    break;
  case 740:
    maximum_value = 0.537489;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.181462, 0.00330738, 0.195896, -0.0548873, 0.108551});
    break;
  case 761:
    maximum_value = 0.373342;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.155269, 0.0176217, 0.152512, -0.0781048, 0.00507746});
    break;
  case 773:
    maximum_value = 0.301145;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.145995, 0.0694747, 0.158823, -0.0379772, -0.0356031});
    break;
  case 793:
    maximum_value = 0.222234;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.116701, 0.0550122, 0.0675169, -0.00402638, -0.012969});
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }
  G4double legendre_cos_theta = 0.;
  for (G4int order = 0, n_order = legendre_coeff.size(); order < n_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
  return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;

} //SigmaZeroPiZeroBC

//______________________________________________________________________________
// calculate diff. cross sec. using of Crystal Ball data (DOI: 10.1103/PhysRevC.80.025204)
// coefficients are used their results
G4bool
SigmaZeroPiZeroCB(const G4double cos_theta, const G4double pk)
{
  std::vector<G4int> pk_list{514, 560, 581, 629, 659, 687, 714, 750};
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

  G4double maximum_value = 0.;  // maximum value of legendre func
  std::vector<G4double> legendre_coeff;
  legendre_coeff.clear();
  switch(pk_min_diff){
  case 514:
    maximum_value = 0.3883;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.1603, -0.0790, 0.0950, -0.0580, 0.0020, 0.0060});
    break;
  case 560:
    maximum_value = 0.3762;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.1598, -0.0874, 0.1120, -0.0350, -0.0020, 0.0160});
    break;
  case 581:
    maximum_value = 0.3958;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.1582, -0.1026, 0.1010, -0.0490, 0.0090, 0.0240});
    break;
  case 629:
    maximum_value = 0.4114;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.1387, -0.1009, 0.1141, -0.0571, 0.0089, 0.0083});
    break;
  case 659:
    maximum_value = 0.3696;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.1397, -0.0897, 0.1202, -0.0410, 0.0076, 0.0286});
    break;
  case 687:
    maximum_value = 0.4025;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.1426, -0.0939, 0.1280, -0.0610, -0.0090, 0.0140});
    break;
  case 714:
    maximum_value = 0.366;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.1437, -0.0547, 0.1530, -0.0460, -0.0044, 0.0270});
    break;
  case 750:
    maximum_value = 0.3409;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.1439, 0.0185, 0.1693, -0.0690, 0.0085, 0.0313});
    break;
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }
  G4double legendre_cos_theta = 0.;
  for (G4int order = 0, n_order = legendre_coeff.size(); order < n_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
  return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
} //SigmaZeroPiZeroCB

//______________________________________________________________________________
// calculate diff. cross sec. using of Bubble Chamber data (https://doi.org/10.1016/0550-3213(70)90461-X)
// coefficients are determined by fitting ourselvs
G4bool
SigmaMinusPiPlus(const G4double cos_theta, const G4double pk)
{
  std::vector<G4int> pk_list{617, 637, 658, 677, 699, 719, 740, 761, 773, 793};
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

  G4double maximum_value = 0.;  // maximum value of legendre func
  std::vector<G4double> legendre_coeff;
  legendre_coeff.clear();
  switch(pk_min_diff){
  case 617:
    maximum_value = 0.621633;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.192099, -0.0252072, 0.33481, -0.00507246, 0.0644446});
    break;
  case 637:
    maximum_value = 0.699038;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.182557, -0.0364784, 0.358166, -0.0623308, 0.0595056});
    break;
  case 658:
    maximum_value = 0.879609;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.200185, -0.0349592, 0.429592, -0.0907064, 0.124166});
    break;
  case 677:
    maximum_value = 0.882071;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.234284, -0.0409681, 0.480279, -0.0783708, 0.0481694});
    break;
  case 699:
    maximum_value = 0.944477;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.240956, -0.0595369, 0.485049, -0.107962, 0.0509729});
    break;
  case 719:
    maximum_value = 1.17816;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.300246, 0.0059601, 0.664913, -0.0754831, 0.14348});
    break;
  case 740:
    maximum_value = 1.00704;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.279039, 0.000700548, 0.560804, -0.112095, 0.0557994});
    break;
  case 761:
    maximum_value = 0.992162;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.278788, 0.0115243, 0.550769, -0.113415, 0.0607135});
    break;
  case 773:
    maximum_value = 0.907868;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.247512, 0.0260744, 0.487977, -0.137496, 0.0609578});
    break;
  case 793:
    maximum_value = 0.60697;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.176799, -0.0219359, 0.30909, -0.107588, -0.00844439});
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }
  G4double legendre_cos_theta = 0.;
  for (G4int order = 0, n_order = legendre_coeff.size(); order < n_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
  return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;

} //SigmaMinusPiPlus

//______________________________________________________________________________
// calculate diff. cross sec. using of Bubble Chamber data (https://doi.org/10.1016/0550-3213(70)90461-X)
// coefficients are determined by fitting ourselvs
G4bool
SigmaPlusPiMinus(const G4double cos_theta, const G4double pk)
{
  std::vector<G4int> pk_list{617, 637, 658, 677, 699, 719, 740, 761, 773, 793};
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

  G4double maximum_value = 0.;  // maximum value of legendre func
  std::vector<G4double> legendre_coeff;
  legendre_coeff.clear();
  switch(pk_min_diff){
  case 617:
    maximum_value = 0.57046;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.297643, -0.231051, 0.0657447, -0.0551197, -0.079333});
    break;
  case 637:
    maximum_value = 0.510569;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.268198, -0.187816, 0.032623, -0.0412552, -0.0193236});
    break;
  case 658:
    maximum_value = 0.536133;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.28283, -0.21724, 0.0506597, -0.0519926, -0.0665898});
    break;
  case 677:
    maximum_value = 0.377612;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.231024, -0.147261, 0.0430769, 0.0216758, -0.0368722});
    break;
  case 699:
    maximum_value = 0.772245;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.237773, -0.256457, 0.174323, -0.0805352, 0.023156});
    break;
  case 719:
    maximum_value = 0.54994;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.205324, -0.226257, 0.118012, -0.04796, -0.0476128});
    break;
  case 740:
    maximum_value = 0.548285;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.185234, -0.122724, 0.153056, -0.139147, -0.0518759});
    break;
  case 761:
    maximum_value = 0.254787;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.170264, -0.000936783, 0.0831127, -0.0751872, -0.121173});
    break;
  case 773:
    maximum_value = 0.338741;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.194971, 0.0734801, 0.0559123, -0.15383, -0.101156});
    break;
  case 793:
    maximum_value = 0.262141;
    legendre_coeff.insert(std::cbegin(legendre_coeff), {0.151723, 0.0707109, 0.0441094, -0.102726, -0.067457});
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }
  G4double legendre_cos_theta = 0.;
  for (G4int order = 0, n_order = legendre_coeff.size(); order < n_order; order++) legendre_cos_theta += legendre_coeff[order]*Kinematics::Legendre(order, cos_theta);
  return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;

} //SigmaPlusPiMinus


} //namespace DiffCrossSection
