#include "CmdLineOpts.h"

CmdLineOpts::CmdLineOpts(int argc, char** argv): argc(argc) {
  std::string opt_flag;
  this->mode = "slam";
  this->out_dir = "./";
  this->group_by = "";
  this->reference = "";
  this->black_list_file = "";

  this->max_mm = 10;
  this->excl_flag = 0;
  this->incl_flag = 0;

  this->num_chunks = 1;
  this->num_threads = 1;
  this->chunk_size = 50000;
  this->load_max = 100;

  this->excl_unpaired = false;

  if (argc < 2) {
    this->help();
    exit(1);
  }

  // Simple loop to parse command line options.
  for (int i = 1; i < argc; i++) {
    opt_flag = argv[i];

    if (argv[i][0] == '-') {
      if (opt_flag == "--mode" || opt_flag == "-m") {
        this->mode = argv[i + 1]; i++;
      } else if (opt_flag == "--group-by" || opt_flag == "-g") {
        this->group_by = argv[i + 1]; i++;
      } else if (opt_flag == "--reference" || opt_flag == "-r") {
        this->reference = argv[i + 1]; i++;
      } else if (opt_flag == "--black-list-file" || opt_flag == "-b") {
        this->black_list_file = argv[i + 1]; i++;
      } else if (opt_flag == "--num-threads" || opt_flag == "-t") {
        this->num_threads = std::stoi(argv[i + 1]); i++;
      } else if (opt_flag == "--num-chunks" || opt_flag == "-c") {
        this->num_chunks = std::stoi(argv[i + 1]); i++;
      } else if (opt_flag == "--chunk-size" || opt_flag == "-s") {
        this->chunk_size = std::stoi(argv[i + 1]); i++;
      } else if (opt_flag == "--excl-flag" || opt_flag == "-F") {
        this->excl_flag = std::stoi(argv[i + 1]); i++;
      } else if (opt_flag == "--incl-flag" || opt_flag == "-f") {
        this->excl_flag = std::stoi(argv[i + 1]); i++;
      } else if (opt_flag == "--max-mm" || opt_flag == "-M") {
        this->max_mm = std::stoi(argv[i + 1]); i++;
      } else if (opt_flag == "--load-max" || opt_flag == "-n") {
        this->load_max = std::stoi(argv[i + 1]); i++;
      } else if (opt_flag == "--excl-unpaired" || opt_flag == "-u") {
        this->excl_unpaired = true;
      }  else {
        if (opt_flag == "--help" || opt_flag == "-h") {
          this->help();
          exit(0);
        } else {
          std::cerr << "Unknown option: " << argv[i] << std::endl;
          this->help();
          exit(1);
        }
      }

    } else {
      if (i == argc - 1) {
        this->bam_file = argv[i];
      } else if (i == argc - 2) {
        this->bam_file = argv[i];
        this->out_dir = argv[i + 1];
      } else {
        this->help();
      }
      break;
    }
  }
}

void CmdLineOpts::help() {
  this->usage();
  std::cerr << std::endl;

  std::cerr << "Help: " << std::endl;
  std::cerr << "  -m, --mode            Running mode [dedup,slam]. Default: slam" << std::endl;
  std::cerr << "  -g, --group_by        Group by column. Default: None" << std::endl;
  std::cerr << "  -r, --reference       Reference file. Default: None" << std::endl;
  std::cerr << "  -b, --black_list_file Black list file. Default: None" << std::endl;
  std::cerr << "  -t, --num-threads     Number of threads. Default: 1" << std::endl;
  std::cerr << "  -c, --num-chunks      Number of chunks. Default: 1" << std::endl;
  std::cerr << "  -s, --chunk-size      Chunk size. Default: 50000" << std::endl;
  std::cerr << "  -F, --excl-flag       Exclude flag. Default: 0" << std::endl;
  std::cerr << "  -f, --incl-flag       Include flag. Default: 0" << std::endl;
  std::cerr << "  -M, --max-mm          Maximum number of mismatches. Default: 0" << std::endl;
  std::cerr << "  -n, --load-max        Load maximum number of records. Default: 100" << std::endl;
  std::cerr << "  -u, --excl-unpaired   Exclude unpaired reads. Default: false" << std::endl;
  std::cerr << "  -h, --help            Print this help message" << std::endl;
}

void CmdLineOpts::usage() {
  std::cerr << "Usage: slamtk [options] <BAM-file> [output-dir]" << std::endl;
}

void CmdLineOpts::print_options() {
  std::cout << "Options: " << std::endl;
  std::cout << "Group reads by: " << this->group_by << std::endl;
}
