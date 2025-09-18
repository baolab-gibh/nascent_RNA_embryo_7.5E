#pragma once
#include <string>
#include <vector>
#include <memory>

#include "CigarString.h"

class Mutations {
private:
  std::vector<std::pair<size_t, std::string>> _pos_mut_vec;

public:
  Mutations();
  Mutations(CigarStringSPtr cigar);
};
