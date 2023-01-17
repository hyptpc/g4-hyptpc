#ifndef SACRunAction_hh
#define SACRunAction_hh

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include "g4root.hh"


#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>

class G4Run;
class SACAnalysisManager;


class SACRunAction : public G4UserRunAction
{
public:
  SACRunAction(SACAnalysisManager *analysisManager = 0);
  virtual ~SACRunAction();

  virtual void BeginOfRunAction(const G4Run *aRun);
  virtual void EndOfRunAction(const G4Run *aRun);

private:
  SACAnalysisManager *anaMan;
};

#endif
