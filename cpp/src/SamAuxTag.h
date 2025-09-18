#pragma once
#include <vector>
#include <string>
#include <memory>
#include <iostream>

#include "Typedefs.h"
#include "htslib/sam.h"

template<class T>
class SamAuxTag {
  std::string name;
  std::string type;
  T value;

public:
  SamAuxTag() = delete;
  SamAuxTag<T>(const std::string&, const std::string&, const T&);
  ~SamAuxTag();

  const std::string& get_name() const;
  const std::string& get_type() const;
  const T& get_value() const;

  const std::string to_string() const;
  int get_val_int() const;
  float get_val_float() const;
  double get_val_double() const;
  uint8_t get_val_uint8() const;
  std::string get_val_string() const;
};

typedef SamAuxTag<std::string> SamAuxTagString;
typedef std::shared_ptr<SamAuxTag<std::string>> SamAuxTagSPtr;

// Force to compile for int, char, float, and double, uint8_t, and std::string.
template class SamAuxTag<int>;
template class SamAuxTag<float>;
template class SamAuxTag<double>;
template class SamAuxTag<std::string>;
