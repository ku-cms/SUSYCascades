#ifndef GitHashTool_H
#define GitHashTool_H

#include <string>

class GitHashTool {
public:
  GitHashTool() = default;

  // Get the current Git commit hash as a string
  std::string GetGitCommitHash() const;

  // Write the Git commit hash to the given ROOT file as a TParameter<std::string>
  void WriteToROOT(const std::string& outputFileName) const;
};

#endif
