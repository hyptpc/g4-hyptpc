#ifndef SACEVENTACTION_HH
#define SACEventAction_hh 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <array>


class SACAnalysisManager;


class SACEventAction : public G4UserEventAction
{
public:
  SACEventAction(SACAnalysisManager *analysisManager=0);
  virtual ~SACEventAction();
private:
  SACEventAction(const SACEventAction &);
  SACEventAction & operator = (const SACEventAction &);

public:
  
  virtual void BeginOfEventAction(const G4Event *anEvent);
  virtual void EndOfEventAction(const G4Event *anEvent);

protected:
  SACAnalysisManager *anaMan;

  
};

#endif

 
  
