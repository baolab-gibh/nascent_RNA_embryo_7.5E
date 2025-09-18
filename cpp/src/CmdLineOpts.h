#pragma once
#include <iostream>
#include <memory>
#include "Typedefs.h"

class CmdLineOpts {
  int argc;

  StrVec argv;
  StrStrMap opts;

  StrVec opts_long{
    "--mode", "--group-by", "--reference", "--black-list-file", "--num-threads", "--num-chunks", "--chunk-size",
    "--excl-flag", "--incl-flag", "--max-mm", "--excl-unpaired", "--load-max"
  };
  StrVec opts_short{"-m", "-g", "-r", "-b", "-t", "-c", "-s"};

public:
  std::string mode;
  std::string out_dir;
  std::string bam_file;
  std::string group_by;
  std::string reference;
  std::string black_list_file;

  size_t max_mm;
  size_t excl_flag;
  size_t incl_flag;

  size_t num_threads;
  size_t num_chunks;
  size_t chunk_size;
  size_t load_max;

  bool excl_unpaired;

public:
  CmdLineOpts(int argc, char** argv);
  ~CmdLineOpts() = default;

  void help();
  void usage();
  void print_options();
};

typedef std::shared_ptr<CmdLineOpts> CmdLineOptsSPtr;
