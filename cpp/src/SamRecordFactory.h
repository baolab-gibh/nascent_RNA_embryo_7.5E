#pragma once
#include <memory>
#include "SamRecord.h"
#include "GenomicRegion.h"
#include "Configs.h"

class SamRecordFactory {
  std::string filename;
  std::vector<SamRecordSPtr> read_pool;
  std::vector<GenomicRegion> genomic_regions;
  unsigned int n_workers;

public:
  SamRecordFactory() = delete;
  SamRecordFactory(const std::string& filename, unsigned int n_workers);
  ~SamRecordFactory() = default;

  void remove_duplicates(ReadsGroupBy g, DedupStrategy d);

  void estimate_t2c_rate();
  void estimate_ntr_ratio();
  void identify_nascent_reads();
  void dump_results(const std::string& out_dir);

private:
  void group_reads_by_umi();
  void group_reads_by_pos();
  void group_reads_by(ReadsGroupBy group_by);
};
