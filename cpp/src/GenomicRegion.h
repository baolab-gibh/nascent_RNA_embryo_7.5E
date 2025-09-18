#pragma once
#include <memory>
#include <string>

class GenomicRegion {
  std::string chrom;
  size_t start;
  size_t end;

public:
  GenomicRegion() = delete;
  GenomicRegion(std::string c, size_t s, size_t e);
  ~GenomicRegion() = default;

  // Getters
  std::string get_chr() const;
  size_t get_start() const;
  size_t get_end() const;

  // Comparisons
  bool operator!=(const GenomicRegion& o) const;
  bool operator==(const GenomicRegion& o) const;
  bool operator< (const GenomicRegion& o) const;
  bool operator<=(const GenomicRegion& o) const;
  bool operator> (const GenomicRegion& o) const;
  bool operator>=(const GenomicRegion& o) const;
  bool operator!=(const std::shared_ptr<GenomicRegion> o) const;
  bool operator==(const std::shared_ptr<GenomicRegion> o) const;
  bool operator< (const std::shared_ptr<GenomicRegion> o) const;
  bool operator<=(const std::shared_ptr<GenomicRegion> o) const;
  bool operator> (const std::shared_ptr<GenomicRegion> o) const;
  bool operator>=(const std::shared_ptr<GenomicRegion> o) const;


  // Operations
  void substract(const GenomicRegion& o);
  void substract(const std::shared_ptr<GenomicRegion> o);
  void merge(const GenomicRegion& o);
  void merge(const std::shared_ptr<GenomicRegion> o);

  //
  GenomicRegion operator-(const GenomicRegion& o);
  GenomicRegion operator+(const GenomicRegion& o);
  GenomicRegion operator-(const std::shared_ptr<GenomicRegion> o);
  GenomicRegion operator+(const std::shared_ptr<GenomicRegion> o);
};
typedef std::shared_ptr<GenomicRegion> GenomicRegionPtr;
