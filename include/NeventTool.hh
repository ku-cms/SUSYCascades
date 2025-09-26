#ifndef NeventTool_h
#define NeventTool_h

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <future>
#include <chrono>
#include <thread>
#include <vector>
#include <unistd.h>     // fork, pipe, read, write, close
#include <sys/wait.h>   // waitpid
#include <signal.h>     // kill
#include <sys/select.h> // select()

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TChain.h"
#include "TError.h"

bool check_dataset_file(std::string dataset_name);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);

class NeventTool {

public:
  NeventTool();
  virtual ~NeventTool();

  void BuildMap(const std::string& rootfile);
  void BuildFilterEffMap(const std::string& rootfile);

  double GetNevent_BKG(const std::string& dataset, const std::string& filetag) const;
  double GetNevent_Cascades(const std::string& dataset, const std::string& filetag) const;
  double GetNevent_SMS(const std::string& dataset, const std::string& filetag, int MP, int MC) const;
  double GetNevent_SMS_code(const std::string& dataset, const std::string& filetag, int MP, int MC, int Code) const;
  double GetNweight_BKG(const std::string& dataset, const std::string& filetag) const;
  double GetNweight_Cascades(const std::string& dataset, const std::string& filetag) const;
  double GetNweight_SMS(const std::string& dataset, const std::string& filetag, int MP, int MC) const;
  double GetNweight_SMS_code(const std::string& dataset, const std::string& filetag, int MP, int MC, int Code) const;
  double GetFilterEff(const std::string& dataset, const std::string& filetag, int lumiblock = -1) const;

  bool DatasetIsFastSim(const std::string& infile);
  int EventsInDAS(const std::string& u_file = "");
  std::string Get_DASdatasetname(const std::string& u_file);
  std::string Get_DASfilename(const std::string& u_file);

private:
  static std::map<std::pair<std::string,std::string>,double> m_Label2Nevent_BKG;
  static std::map<std::pair<std::string,std::string>,double> InitMap_Nevent_BKG();
  static std::map<std::pair<std::string,std::string>,double> m_Label2Nweight_BKG;
  static std::map<std::pair<std::string,std::string>,double> InitMap_Nweight_BKG();
  static std::map<std::pair<std::string,std::string>,double> m_Label2Nevent_Cascades;
  static std::map<std::pair<std::string,std::string>,double> InitMap_Nevent_Cascades();
  static std::map<std::pair<std::string,std::string>,double> m_Label2Nweight_Cascades;
  static std::map<std::pair<std::string,std::string>,double> InitMap_Nweight_Cascades();
  static std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > m_Label2Nevent_SMS;
  static std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > InitMap_Nevent_SMS();
  static std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > m_Label2Nweight_SMS;
  static std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > InitMap_Nweight_SMS();
  static std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> m_Label2Nevent_SMS_code;
  static std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> InitMap_Nevent_SMS_code();
  static std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> m_Label2Nweight_SMS_code;
  static std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> InitMap_Nweight_SMS_code();

  static std::map<std::string,std::map<int,double> > InitMap_FilterEff();
  static std::map<std::string,std::map<int,double> > m_Label2FilterEff;

  std::string m_RootFile;
  TFile* m_File;
  TTree* m_Tree;
  std::string m_RootFile_FE;
  TFile* m_File_FE;
  TTree* m_Tree_FE;

  int m_lumiblock;
  double m_efficiency;
  
  double m_Nevent;
  int m_NDAS;
  double m_Nweight;
  std::string* m_dataset;
  std::string* m_filetag;
  int m_MP;
  int m_MC;
  int m_Code;
  TBranch* b_m_lumiblock;  
  TBranch* b_m_efficiency; 
  TBranch* b_m_NDAS;  
  TBranch* b_m_Nevent;  
  TBranch* b_m_Nweight; 
  TBranch* b_m_dataset;
  TBranch* b_m_filetag;
  TBranch* b_m_MP;
  TBranch* b_m_MC;
  TBranch* b_m_Code;

  void Initialize_BKG(const std::string& dataset, const std::string& filetag) const;
  void Initialize_Cascades(const std::string& dataset, const std::string& filetag) const;
  void Initialize_SMS(const std::string& dataset, const std::string& filetag) const;
  void Initialize_SMS_code(const std::string& dataset, const std::string& filetag) const;
  void Initialize_FE(const std::string& dataset, const std::string& filetag) const;
  
};

// Helper: build a new TChain and add files, optionally replacing a src prefix with dstPrefix.
inline TChain* buildChain(const std::vector<std::string>& filenames,
                          const std::string& srcPrefix,
                          const std::string& dstPrefix,
                          const char* treeName,
                          bool useTree)
{
    TChain* tmp = nullptr;
    if (useTree) tmp = new TChain(treeName);
    else tmp = new TChain("Events");

    for (const auto& f : filenames) {
        std::string fname = f;
        if (!dstPrefix.empty() && !srcPrefix.empty()) {
            auto pos = fname.find(srcPrefix);
            if (pos != std::string::npos) {
                fname.replace(pos, srcPrefix.size(), dstPrefix);
            }
        }
        tmp->Add(fname.c_str());
        std::cout << "   Added file " << fname << std::endl;
    }
    return tmp;
}

// Helper: run chain->GetEntries() in a subprocess with timeout (seconds).
// Returns -1 on timeout or error, otherwise the entry count.
inline Long64_t tryGetEntriesSubproc(TChain* chain, int timeoutSeconds = 60)
{
    int prevErrLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;

    int pipefd[2];
    if (pipe(pipefd) != 0) {
        perror("pipe");
        gErrorIgnoreLevel = prevErrLevel;
        return -1;
    }

    pid_t pid = fork();
    if (pid < 0) {
        perror("fork");
        gErrorIgnoreLevel = prevErrLevel;
        return -1;
    }

    if (pid == 0) {
        // Child process
        close(pipefd[0]); // close read end
        Long64_t entries = chain->GetEntries();
        if (write(pipefd[1], &entries, sizeof(entries)) < 0) {
            // ignore write error
        }
        close(pipefd[1]);
        _exit(0);
    }

    // Parent process
    close(pipefd[1]); // close write end
    Long64_t result = -1;

    fd_set set;
    FD_ZERO(&set);
    FD_SET(pipefd[0], &set);

    struct timeval tv;
    tv.tv_sec = timeoutSeconds;
    tv.tv_usec = 0;

    int rv = select(pipefd[0] + 1, &set, nullptr, nullptr, &tv);
    if (rv > 0) {
        if (read(pipefd[0], &result, sizeof(result)) < 0) {
            perror("read");
            result = -1;
        }
        int status = 0;
        waitpid(pid, &status, 0);
    } else {
        // timeout or error
        kill(pid, SIGKILL);
        int status = 0;
        waitpid(pid, &status, 0);
        result = -1;
    }

    close(pipefd[0]);
    gErrorIgnoreLevel = prevErrLevel;
    return result;
}

#endif
