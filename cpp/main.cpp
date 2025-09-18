#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <memory>
#include <algorithm>
#include <utility>

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/hts_endian.h"

#include "SamRecord.h"
#include "SamAuxTag.h"
#include "BamReader.h"
#include "CmdLineOpts.h"
#include "Typedefs.h"

// typedef std::string Str;
typedef std::shared_ptr<std::string> StrSPtr;
typedef struct kstring_t KString;

void run(CmdLineOptsSPtr opts) {
  BamReader::BamReaderSPtr bam_file(new BamReader(opts->bam_file));
  for (size_t x = 0; x < opts->load_max; x++) {
    auto rec = bam_file->next();

    if (rec->skipRead(opts->excl_flag, opts->incl_flag, opts->max_mm, opts->excl_unpaired)) {
      continue;
    }

    // if (rec->hasT2C(100, 0, 0)) { std::cout << rec << std::endl; }

    if (rec == NULL) { break; }
  }
}

int main(int argc, char** argv) {
  CmdLineOptsSPtr opts(new CmdLineOpts(argc, argv));
  run(opts);

  return 0;
}
