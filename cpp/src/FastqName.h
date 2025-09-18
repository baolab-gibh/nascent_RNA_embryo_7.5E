#pragma once
#include <string>
#include <memory>

class FastqName {
  std::string name;

public:
  FastqName(const std::string&);
  ~FastqName();

  const std::string& guess_version();
  const std::string& to_string() const { return this->name; }
};
