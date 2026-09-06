// -*- C++ -*-

#include "ConfMan.hh"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <iterator>
#include <libgen.h>
#include <sstream>
#include <string>
#include <vector>

#include <CLHEP/Units/SystemOfUnits.h>

#include <TFile.h>
#include <TNamed.h>

#include "BeamMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "HistMan.hh"
#include "JamMan.hh"
#include "IncMan.hh"
#include "DiffCrossSectionMan.hh"

//_____________________________________________________________________________
ConfMan::ConfMan()
  : m_conf_key("CONF"),
    m_conf_dir(),
    m_is_ready(false),
    m_file(),
    m_string(),
    m_double(),
    m_int(),
    m_bool()
{
}

//_____________________________________________________________________________
ConfMan::~ConfMan()
{
}

//_____________________________________________________________________________
G4bool
ConfMan::Contains(const G4String& key) const
{
  return m_string.find(key) != m_string.end();
}

//_____________________________________________________________________________
std::set<G4int>
ConfMan::GetIntList(const G4String& key) const
{
  std::set<G4int> out;
  auto it = m_string.find(key);
  if (it == m_string.end()) return out;

  std::stringstream ss(it->second.data());
  G4int v = 0;
  while (ss >> v) {
    out.insert(v);
  }
  return out;
}

//_____________________________________________________________________________
std::set<G4int>
ConfMan::GetOrDefaultIntList(const G4String& key,
                             const std::set<G4int>& default_value) const
{
  return Contains(key) ? GetIntList(key) : default_value;
}

//_____________________________________________________________________________
G4bool
ConfMan::Initialize()
{
  if(m_is_ready){
    G4cerr << FUNC_NAME
	   << " already initialied" << G4endl;
    return false;
  }

  std::ifstream ifs(m_file[m_conf_key]);
  if(!ifs.is_open()){
    G4cerr << FUNC_NAME
	   << " cannot open file : " << m_file[m_conf_key] << G4endl;
    return false;
  }

  G4cout << FUNC_NAME << G4endl
	 << " open file : " << m_file[m_conf_key] << G4endl;
  m_string[m_conf_key] = m_file[m_conf_key];

  m_conf_dir = ::dirname(const_cast<char*>(m_file[m_conf_key].data()));
  m_conf_buf.clear();
  m_conf_buf += "\n";

  G4String line;
  while(ifs.good() && std::getline(ifs, line)){
    m_conf_buf += line + "\n";
    if (line[0]=='#')
      continue;
    std::istringstream iss(line);
    G4String key;
    iss >> key;
    if (key.empty())
      continue;
    if (key.back() == ':')
      key.pop_back();

    // Rest of the line is the value (allows space/tab-separated lists).
    std::string rest;
    std::getline(iss >> std::ws, rest);
    // Drop inline comment.
    const auto hash = rest.find('#');
    if (hash != std::string::npos)
      rest.erase(hash);
    // Trim trailing whitespace.
    while (!rest.empty() && std::isspace(static_cast<unsigned char>(rest.back())))
      rest.pop_back();
    G4String val = rest;
    if (val.empty())
      continue;

    G4cout << " key = "   << std::setw(20) << std::left << key
	   << " value = " << std::setw(30) << std::left << val
	   << G4endl;
    m_file[key]   = FilePath(val);
    m_string[key] = val;
    m_double[key] = std::strtod(val, nullptr);
    m_int[key]    = std::strtol(val, nullptr, 10);
    m_bool[key]   = static_cast<G4bool>(std::strtol(val, nullptr, 10));
  }

  if (!InitializeParameterFiles() || !InitializeHistograms())
    return false;

  // if(gUser.IsReady())
  //   gUser.Print();

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
ConfMan::Initialize(const G4String& file_name)
{
  m_file[m_conf_key] = file_name;
  return Initialize();
}

//_____________________________________________________________________________
G4bool
ConfMan::Initialize(const std::vector<G4String>& arg)
{
  m_file[m_conf_key] = arg[kConfFile];
  m_string["ROOT"] = arg[kOutFile];
  return Initialize();
}

//_____________________________________________________________________________
G4bool
ConfMan::InitializeHistograms()
{
  return true;
}

//_____________________________________________________________________________
G4bool
ConfMan::InitializeParameterFiles()
{
  return (true
          && InitializeParameter<DCGeomMan>("DCGEO")
          && InitializeParameter<BeamMan>("BEAM")
          && InitializeParameter<DetSizeMan>("DSIZE")
          && InitializeParameter<HistMan>("HIST")
          && InitializeParameter<JamMan>("JAM")
          && InitializeParameter<IncMan>("INC")
	  && InitializeParameter<DiffCrossSectionMan>("DCS")
          );
}

//_____________________________________________________________________________
// G4bool
// ConfMan::Finalize()
// {
//   return FinalizeProcess();
// }

//_____________________________________________________________________________
G4String
ConfMan::FilePath(const G4String& src) const
{
  std::ifstream tmp(src);
  if (tmp.good())
    return src;
  else
    return m_conf_dir + "/" + src;
}
