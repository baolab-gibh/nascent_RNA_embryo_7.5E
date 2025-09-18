#include "GenomicRegion.h"


inline bool GenomicRegion::operator==(const GenomicRegion& o) const {
  return this->chrom == o.chrom && this->start == o.start && this->end == o.end;
}
inline bool GenomicRegion::operator!=(const GenomicRegion& o) const {
  return !(*this == o);
}
inline bool GenomicRegion::operator< (const GenomicRegion& o) const {
  return this->chrom == o.chrom ? this->start < o.start : this->chrom < o.chrom;
}
inline bool GenomicRegion::operator<=(const GenomicRegion& o) const {
  return *this < o || *this == o;
}
inline bool GenomicRegion::operator> (const GenomicRegion& o) const {
  return this->chrom == o.chrom ? this->start > o.start : this->chrom > o.chrom;
}
inline bool GenomicRegion::operator>=(const GenomicRegion& o) const {
  return *this > o || *this == o;
}

inline bool GenomicRegion::operator==(const GenomicRegionPtr o) const {
  return this->chrom == o->chrom && this->start == o->start && this->end == o->end;
}
inline bool GenomicRegion::operator!=(const GenomicRegionPtr o) const {
  return !(*this == o);
}
inline bool GenomicRegion::operator< (const GenomicRegionPtr o) const {
  return this->chrom == o->chrom ? this->start < o->start : this->chrom < o->chrom;
}
inline bool GenomicRegion::operator<=(const GenomicRegionPtr o) const {
  return *this < o || *this == o;
}
inline bool GenomicRegion::operator> (const GenomicRegionPtr o) const {
  return this->chrom == o->chrom ? this->start > o->start : this->chrom > o->chrom;
}
inline bool GenomicRegion::operator>=(const GenomicRegionPtr o) const {
  return *this > o || *this == o;
}

// Getters
std::string GenomicRegion::get_chr() const { return this->chrom; }
size_t GenomicRegion::get_start() const { return this->start; }
size_t GenomicRegion::get_end() const { return this->end; }

// Operations
void GenomicRegion::substract(const GenomicRegion& o) {
  if (this->chrom == o.chrom) {
    auto is_overlapped = this->start <= o.end && this->end >= o.start;
  } else {
    return;
  }
}

void GenomicRegion::substract(const GenomicRegionPtr o) {
  if (this->chrom == o->chrom) {
    auto is_overlapped = this->start <= o->end && this->end >= o->start;
  } else {
    return;
  }
}

void GenomicRegion::merge(const GenomicRegion& o) {
  if (this->chrom == o.chrom) {
    auto is_overlapped = this->start <= o.end && this->end >= o.start;
  } else {
    return;
  }
}

void GenomicRegion::merge(const GenomicRegionPtr o) {
  if (this->chrom == o->chrom) {
    auto is_overlapped = this->start <= o->end && this->end >= o->start;
  } else {
    return;
  }
}
