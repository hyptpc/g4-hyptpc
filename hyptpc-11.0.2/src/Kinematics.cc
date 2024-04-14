// -*- C++ -*-

#include "Kinematics.hh"

#include <Randomize.hh>

#include "AnaManager.hh"
#include "FuncName.hh"

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
}


namespace Kinematics
{

const G4double b = 1.62 * CLHEP::GeV * CLHEP::fermi / CLHEP::hbarc;


//_____________________________________________________________________________
G4ThreeVector
HarmonicFermiMomentum(G4int type)
{
  /*      THIS ROUTINE GENERATES FERMI MOMENTUM KF(3) BY HARMONIC */
  /*      OSCILATOR MODEL FOR 1S AND 1P */
  /*      INPUT; L=0 OR 1 */
  /*      OUTPUT; KF(3) */
  ///what is unit?? mev? gev?
  ////in the g3LEPS, b=1.64 (fm)
  ////in the g3LEPS, b=1.73 hicks
  ////in the g3LEPS, b=1.93 fitting
  ////hbar = 197.327053 MeV*fm --. 0.197 --> GeV*fm

  G4double x = 0.;
  switch(type){
  case 0: {
    static const G4double ymax = std::exp(-1);
    while(true){
      x = G4RandFlat::shoot();
      G4double y = x * x * std::exp(-b * b * x * x);
      G4double yy = G4RandFlat::shoot() * ymax;
      // G4cout << "iterate " << x << " " << b << " "
      // 	     << y << " " << yy << G4endl;
      if (yy < y) break;
    }
    break;
  }
  case 1: {
    static const G4double ymax = std::exp(-2) * 4 / (b * b);
    while(true){
      x = G4RandFlat::shoot();
      G4double y = b * b * x * x * x * x * std::exp(-b * b * x * x);
      G4double yy = G4RandFlat::shoot() * ymax;
      // G4cout << "iterate " << x << " " << b << " "
      // 	     << y << " " << yy << G4endl;
      if (yy < y) break;
    }
    break;
  }
  default:
    G4Exception(FUNC_NAME, " ", RunMustBeAborted, "");
    break;
  }

  G4double theta = std::acos(G4RandFlat::shoot(-1., 1.));
  G4double phi   = G4RandFlat::shoot(-CLHEP::pi, CLHEP::pi);

  return G4ThreeVector(x * std::sin(theta) * std::cos(phi),
			x * std::sin(theta) * std::sin(phi),
			x * std::cos(theta));
}

//______________________________________________________________________________
G4int
HarmonicFermiMomentumDeuteron(G4double* Kf)
{
  //deut_pro : 0: deuteron
  //deut_pro : 1: proton
  //// for deuteron E27
  /////// Copy of H. Takahashi-san's code
  /////// revised to geant4 by S.Hwang
  G4double ymax;
  G4double x, y, yy, theta, phi;
  /* THIS ROUTINE GENERATES FERMI MOMENTUM KF(3) BY HARMONIC */
  /* OSCILATOR MODEL FOR 1S AND 1P */
  /* INPUT; L=0 OR 1 */
  /* OUTPUT; KF(3) */
  ///what is unit?? mev? gev?
  ////in the g3LEPS, b=1.64 (fm)
  ////in the g3LEPS, b=1.73 hicks
  ////in the g3LEPS, b=1.93 fitting
  // G4double b = 1.62;
  // G4double b = 1.93;
  // G4double b = 2.1;

  ymax = 20.;
  do {
    x = G4RandFlat::shoot();
    y = (x*x) * (13.* std::exp(-8.*x)+30000.*std::exp(-37.*x));
    yy = G4RandFlat::shoot() * ymax;
  } while (yy > y);

  theta = std::acos(G4RandFlat::shoot(-1., 1.));
  phi   = (G4RandFlat::shoot(-1., 1.))*CLHEP::pi;
  // IsotropicAngle(&theta, &phi);
  Kf[0] = x * std::sin(theta) * std::cos(phi);
  Kf[1] = x * std::sin(theta) * std::sin(phi);
  Kf[2] = x * std::cos(theta);
  // G4cout<<x<<":"<<Kf[0]<<":"<<Kf[1]<<":"<<Kf[2]<<":"<<G4endl;
  return 0;
}

//______________________________________________________________________________
G4double
Legendre(G4int order, G4double x)
{
  switch(order){
  case 0:
    return 1.;
  case 1:
    return x;
  case 2:
    return 1./2.*(3.*x*x - 1.);
  case 3:
    return 1./2.*(5.*x*x*x - 3.*x);
  case 4:
    return 1./8.*(35.*x*x*x*x - 30.*x*x + 3.);
  case 5:
    return 1./8.*(63.*x*x*x*x*x - 70.*x*x*x + 15.*x);
  case 6:
    return 1./16.*(231.*x*x*x*x*x*x - 315.*x*x*x*x + 105.*x*x - 5.);
  case 7:
    return 1./16.*(429.*x*x*x*x*x*x*x -
		    693.*x*x*x*x*x +
		    315.*x*x*x -
		     35.*x);
  case 8:
    return 1./128.*( 6435.*x*x*x*x*x*x*x*x -
		     12012.*x*x*x*x*x*x +
		      6930.*x*x*x*x -
		      1260.*x*x +
		        35.);
  case 9:
    return 1./128.*(12155.*x*x*x*x*x*x*x*x*x -
		     25740.*x*x*x*x*x*x*x +
		     18018.*x*x*x*x*x -
		      4620.*x*x*x +
		       315.*x);
  case 10:
    return 1./256.*( 46189.*x*x*x*x*x*x*x*x*x*x -
		     109395.*x*x*x*x*x*x*x*x +
		      90090.*x*x*x*x*x*x -
		      30030.*x*x*x*x +
		       3465.*x*x -
		         63.);
  case 11:
    return 1./256.*( 88179.*x*x*x*x*x*x*x*x*x*x*x -
		     230945.*x*x*x*x*x*x*x*x*x +
		     218790.*x*x*x*x*x*x*x -
		      90090.*x*x*x*x*x +
		      15015.*x*x*x -
		        693.*x);
  case 12:
    return 1./1024.*( 676039.*x*x*x*x*x*x*x*x*x*x*x*x -
		      1939938.*x*x*x*x*x*x*x*x*x*x +
		      2078505.*x*x*x*x*x*x*x*x -
		      1021020.*x*x*x*x*x*x +
		       225225.*x*x*x*x -
		        18018.*x*x +
		          231.);
  default:
    G4cout << "#E " << FUNC_NAME << " invalid order : " << order << G4endl;
    return 0.;
  }
}

//______________________________________________________________________________
// calculation the intersection point of circle and line
//   line: x-x1 = u*(z-z1) -> x = u*z + w, (w = x1-u*z1)
//   circ: x^2 + (z-z0)^2 = r^2
//   beam hit pos_in:  (x1, y1, z1)
//   target center:    (x0, y0, z0)
//   beam hit pos_out: (x2, y2, z2)
G4double
EffectiveThickness(const G4ThreeVector pos, const G4ThreeVector mom, const G4ThreeVector target_pos, const G4ThreeVector target_size)
{
  const G4double u = mom.getX()/mom.getZ();
  const G4double v = mom.getY()/mom.getZ();
  G4double x0 = target_pos.getX();
  G4double z0 = target_pos.getZ();
  G4double x1 = pos.getX();
  G4double y1 = pos.getY();
  G4double z1 = pos.getZ();
  G4double  w = x1 - u*z1;
  G4double  r = target_size.getY()/2; 
  G4double sqrt_term = std::sqrt( r*r*(u*u+1) - (u*z0+w-x0)*(u*z0+w-x0) );
  G4double z2 = ( u*(x0-w) + z0 + sqrt_term )/(u*u+1);
  G4double x2 = u*z2 + w;
  G4double y2 = v*(z2-z1) + y1;
  return std::sqrt( (x2-x1)*(x2-x1) + (z2-z1)*(z2-z1) );
  // return std::sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
}

//______________________________________________________________________________
// determine vertex position randomly
// cal. method is same as EffectiveThickness
G4ThreeVector
RandomVertex(const G4ThreeVector pos, const G4ThreeVector mom, const G4ThreeVector target_pos, const G4ThreeVector target_size)
{
  const G4double u = mom.getX()/mom.getZ();
  const G4double v = mom.getY()/mom.getZ();
  G4double x0 = target_pos.getX();
  G4double y0 = target_pos.getY();
  G4double z0 = target_pos.getZ();
  G4double x1 = pos.getX();
  G4double y1 = pos.getY();
  G4double z1 = pos.getZ();
  G4double  w = x1 - u*z1;
  G4double  r = target_size.getY()/2; 
  G4double  h = target_size.getZ()/2; 
  G4double sqrt_term = std::sqrt( r*r*(u*u+1) - (u*z0+w-x0)*(u*z0+w-x0) );
  G4double z2_minus  = ( u*(x0-w) + z0 - sqrt_term )/(u*u+1);
  G4double z2_plus   = ( u*(x0-w) + z0 + sqrt_term )/(u*u+1);

  G4double rand_z = G4RandFlat::shoot(z2_minus, z2_plus);
  G4double rand_x = u*rand_z + w;
  G4double rand_y = v*(rand_z-z1)+y1;
  G4ThreeVector vertex(rand_x, rand_y, rand_z);
  while ( std::sqrt((rand_x-x0)*(rand_x-x0) + (rand_z-z0)*(rand_z-z0)) > r || std::abs(rand_y-y0) > h ) {
    rand_z = G4RandFlat::shoot(z2_minus, z2_plus);
    rand_x = u*rand_z + w;
    rand_y = v*(rand_z-z1)+y1;
    vertex.set(rand_x, rand_y, rand_z);
  }
  return vertex;
}

//______________________________________________________________________________
// calculate diff. cross sec. using CB data
G4bool
CrystalBallLegendre(const G4double cos_theta, const G4double pk)
{
  G4int pk_list[15] = {724, 726, 728, 730, 732, 734, 738, 742, 746, 750, 754, 758, 762, 766, 770};
  G4int pk_min_diff = pk_list[0];
  G4double min_diff = 1000;
  for (G4int i = 0; i < 15; i++) {
    G4double diff = std::abs( pk-pk_list[i] );
    if (diff < min_diff) {
      min_diff    = diff;
      pk_min_diff = pk_list[i];
    }
  }
  gAnaMan.SetDiffCrossMom(pk_min_diff);

  G4int legendre_order = 3;
  switch(pk_min_diff){
  case 724:
    {
      G4double legendre_coeff[legendre_order] = {0.0162793, 0.00440176, 0.00677799};
      G4double maximum_value = 0.0274591;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 726:
    {
      G4double legendre_coeff[legendre_order] = {0.0297064, -0.000898291, 0.0138006};
      G4double maximum_value =  0.0444053;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 728:
    {
      G4double legendre_coeff[legendre_order] = {0.0445592, -0.00801745, 0.0139182};
      G4double maximum_value = 0.0664949;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 730:
    {
      G4double legendre_coeff[legendre_order] = {0.0657561, 0.0202863, 0.0526398};
      G4double maximum_value = 0.138682;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 732:
    {
      G4double legendre_coeff[legendre_order] = {0.0688638, -0.00504829, 0.0363991};
      G4double maximum_value = 0.110311;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 734:
    {
      G4double legendre_coeff[legendre_order] = {0.0902174, 0.0288109, 0.104857};
      G4double maximum_value = 0.223885;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
     }
  case 738:
    {
      G4double legendre_coeff[legendre_order] = {0.10436, 0.0235637, 0.0526217};
      G4double maximum_value = 0.180545;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 742:
    {
      G4double legendre_coeff[legendre_order] = {0.111125, -0.0107509, 0.0518198};
      G4double maximum_value = 0.173696;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 746:
    {
      G4double legendre_coeff[legendre_order] = {0.113412, -0.0159145, 0.0277801};
      G4double maximum_value = 0.157107;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 750:
    {
      G4double legendre_coeff[legendre_order] = {0.0944732, -0.000168073, 0.0425309};
      G4double maximum_value = 0.137172;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 754:
    {
      G4double legendre_coeff[legendre_order] = {0.101783, -0.0137408, 0.0162012};
      G4double maximum_value = 0.131725;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 758:
    {
      G4double legendre_coeff[legendre_order] = {0.0942948, -0.0109754, 0.0459009};
      G4double maximum_value = 0.151171;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 762:
    {
      G4double legendre_coeff[legendre_order] = {0.0723755, -0.044794, -0.00302102};
      G4double maximum_value = 0.114148;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 766:
    {
      G4double legendre_coeff[legendre_order] = {0.0787211, -0.0477353, 0.0156414};
      G4double maximum_value = 0.142098;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  case 770:
    {
      G4double legendre_coeff[legendre_order] = {0.0586327, 0.00726473, 0.0193116};
      G4double maximum_value = 0.0852091;  // maximum value of legendre func
      G4double legendre_cos_theta = 0.;
      for (G4int order = 0; order < legendre_order; order++) legendre_cos_theta += legendre_coeff[order]*Legendre(order, cos_theta);
      return G4RandFlat::shoot(0.0, maximum_value) <= legendre_cos_theta;
    }
  default:
    G4cout << "#E " << FUNC_NAME << " invalid pk_min_diff : " << pk_min_diff << G4endl;
    return false;
  }
}

}
