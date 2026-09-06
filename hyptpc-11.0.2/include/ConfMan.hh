// -*- C++ -*-

#ifndef CONF_MAN_HH
#define CONF_MAN_HH

#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <globals.hh>

enum EArgv { kProcess, kConfFile, kOutFile, kG4Macro, kArgc };

//_____________________________________________________________________________
class ConfMan
{
public:
  static G4String  ClassName();
  static ConfMan& GetInstance();
  ~ConfMan();

private:
  ConfMan();
  ConfMan(const ConfMan&);
  ConfMan& operator=(const ConfMan&);

private:
  typedef std::map<G4String, G4String> StrList;
  typedef std::map<G4String, G4double> DoubleList;
  typedef std::map<G4String, G4int>    IntList;
  typedef std::map<G4String, G4bool>   BoolList;
  G4String     m_conf_key;
  G4String     m_conf_dir;
  G4String     m_conf_buf;
  G4bool       m_is_ready;
  StrList      m_file;
  StrList      m_string;
  DoubleList   m_double;
  IntList      m_int;
  BoolList     m_bool;

public:
  // G4bool    Finalize();
  // G4bool    FinalizeProcess();
  G4bool    Contains(const G4String& key) const;
  template <typename T>
  static const T& Get(const G4String& key);
  template <typename T>
  static T GetOrDefault(const G4String& key, const T& default_value);
  // Space/tab-separated integers from the conf value string.
  std::set<G4int> GetIntList(const G4String& key) const;
  std::set<G4int> GetOrDefaultIntList(const G4String& key,
                                      const std::set<G4int>& default_value) const;
  const G4String& ConfBuf() const { return m_conf_buf; }
  G4bool    Initialize();
  G4bool    Initialize(const G4String& file_name);
  G4bool    Initialize(const std::vector<G4String>& arg);
  G4bool    InitializeHistograms();
  G4bool    InitializeParameterFiles();
  template <typename T>
  G4bool    InitializeParameter();
  template <typename T>
  G4bool    InitializeParameter(const G4String& key);
  template <typename T>
  G4bool    InitializeParameter(const G4String& key1,
                                const G4String& key2);
  G4bool    IsReady() const { return m_is_ready; }

private:
  G4String FilePath(const G4String& src) const;
  G4bool   ShowResult(G4bool status, const G4String& name) const;
};

//_____________________________________________________________________________
inline G4String
ConfMan::ClassName()
{
  static const G4String s_name("ConfMan");
  return s_name;
}

//_____________________________________________________________________________
inline ConfMan&
ConfMan::GetInstance()
{
  static ConfMan s_instance;
  return s_instance;
}

//_____________________________________________________________________________
template <>
inline const G4String&
ConfMan::Get<G4String>(const G4String& key)
{
  return GetInstance().m_string[key];
}

//_____________________________________________________________________________
template <>
inline const G4double&
ConfMan::Get<G4double>(const G4String& key)
{
  return GetInstance().m_double[key];
}

//_____________________________________________________________________________
template <>
inline const G4int&
ConfMan::Get<G4int>(const G4String& key)
{
  return GetInstance().m_int[key];
}

//_____________________________________________________________________________
template <>
inline const G4bool&
ConfMan::Get<G4bool>(const G4String& key)
{
  return GetInstance().m_bool[key];
}

//_____________________________________________________________________________
template <>
inline G4String
ConfMan::GetOrDefault<G4String>(const G4String& key, const G4String& default_value)
{
  return GetInstance().Contains(key) ? GetInstance().m_string.at(key) : default_value;
}

//_____________________________________________________________________________
template <>
inline G4double
ConfMan::GetOrDefault<G4double>(const G4String& key, const G4double& default_value)
{
  return GetInstance().Contains(key) ? GetInstance().m_double.at(key) : default_value;
}

//_____________________________________________________________________________
template <>
inline G4int
ConfMan::GetOrDefault<G4int>(const G4String& key, const G4int& default_value)
{
  return GetInstance().Contains(key) ? GetInstance().m_int.at(key) : default_value;
}

//_____________________________________________________________________________
template <>
inline G4bool
ConfMan::GetOrDefault<G4bool>(const G4String& key, const G4bool& default_value)
{
  return GetInstance().Contains(key) ? GetInstance().m_bool.at(key) : default_value;
}

//_____________________________________________________________________________
inline G4bool
ConfMan::ShowResult(G4bool status, const G4String& name) const
{
  if(status)
    std::cout << std::setw(20) << std::left
	      << " ["+name+"]"
	      << "-> Initialized" << std::endl;
  else
    std::cout << std::setw(20) << std::left
	      << " ["+name+"]"
	      << "-> Failed" << std::endl;
  return status;
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter()
{
  return
    ShowResult(T::GetInstance().Initialize(),
               T::GetInstance().ClassName());
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter(const G4String& key)
{
  return
    ShowResult(T::GetInstance().Initialize(m_file[key]),
               T::GetInstance().ClassName());
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter(const G4String& key1,
                             const G4String& key2)
{
  return
    ShowResult(T::GetInstance().Initialize(m_file[key1],
                                           m_file[key2]),
               T::GetInstance().ClassName());
}

#endif
