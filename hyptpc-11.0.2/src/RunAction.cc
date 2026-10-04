// -*- C++ -*-

#include "RunAction.hh"

#include <fstream>
#include <random>

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4StateManager.hh>
#include <G4Timer.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>
#include <Randomize.hh>

#include "AnaManager.hh"
#include "ConfMan.hh"
#include "FuncName.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();
auto& gAnaMan = AnaManager::GetInstance();
G4Timer timer;
}

//_____________________________________________________________________________
RunAction::RunAction()
  : G4UserRunAction()
{
}

//_____________________________________________________________________________
RunAction::~RunAction()
{
}

//_____________________________________________________________________________
void
RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << FUNC_NAME << G4endl
	 << "   Run# = " << aRun->GetRunID() << G4endl;
  gAnaMan.BeginOfRunAction(aRun->GetRunID());

  G4long seed;
  if (gConf.Contains("RandomSeed")) {
    seed = gConf.Get<G4int>("RandomSeed");
  } else {
    std::random_device rd;
    seed = static_cast<G4long>(rd());
    if (seed == 0) seed = 1;
  }
  seed += static_cast<G4long>(aRun->GetRunID());
  G4Random::setTheSeed(seed);
#ifdef DEBUG
  G4cout << "   Initial Seed = " << G4Random::getTheSeed() << G4endl;
  G4Random::showEngineStatus();
#endif
  timer.Start();
}

//_____________________________________________________________________________
void
RunAction::EndOfRunAction(const G4Run* aRun)
{
  timer.Stop();
  gAnaMan.EndOfRunAction();
  G4cout << FUNC_NAME << G4endl
	 << "   Process end  = " << timer.GetClockTime()
	 << "   Event number = " << aRun->GetNumberOfEvent() << G4endl
	 << "   Elapsed time = " << timer << G4endl << G4endl;
}
