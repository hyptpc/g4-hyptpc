#ifndef SACAnalysisManager_h
#define SACAnalysisManager_h

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <cmath>
#include <vector>

class G4Run;
class G4Event;
class TFile;
class TTree;
class EvtVector4R;

const int num_evtgen = 100000;
const int num_aerohit = 100000;
const int num_pmt15hit = 100000;


class SACAnalysisManager
{
public:
  void BeginOfRun(const G4Run*);
  void EndOfRun(const G4Run*);
  void BeginOfEvent(const G4Event*);
  void EndOfEvent(const G4Event*);

private:
  G4String outfile;
  TFile *hfile;
  TTree *tree;

public:
  SACAnalysisManager(const G4String &histname);
  virtual ~SACAnalysisManager();
  SACAnalysisManager(const SACAnalysisManager&);
  SACAnalysisManager& operator=(const SACAnalysisManager&);

  void SetEvtGen(int j, int partnum);
  void Terminate(void) const;
  void SaveFile(void) const;

private:
  G4bool fActive_;
  int event;

  int nEvt;
  int evtid[num_evtgen];
  //int evtpid[num_evtgen];
  int evtpid;
  double evtposx;
  double evtposy;
  double evtposz;
  int evtnumce;
  double evtwavelength[num_evtgen];
  


  int nhAero;
  double aeroangle[num_aerohit];
  int aeropid[num_aerohit];
  double aerotime[num_aerohit];
  double aeroposx[num_aerohit];
  double aeroposy[num_aerohit];
  double aeroposz[num_aerohit];

  int nhPmt15;
  int pmt15pid[num_pmt15hit];
  double pmt15time[num_pmt15hit];
  double pmt15posx[num_pmt15hit];
  double pmt15posy[num_pmt15hit];
  double pmt15posz[num_pmt15hit];
  double pmt15wavelength[num_pmt15hit];
  int pmt15num[num_pmt15hit];


};

#endif
