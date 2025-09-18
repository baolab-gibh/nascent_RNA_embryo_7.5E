#include <string>
#include <vector>
#include <array>
#include <map>
#include "htslib/bgzf.h"

typedef std::map<std::string, size_t> StrSizeTMap;
typedef std::map<std::string, std::string> StrStrMap;

typedef std::pair<std::string, size_t> StrSizeTPair;
typedef std::pair<std::string, std::string> StrStrPair;

typedef std::vector<std::string> StrVec;
typedef std::vector<std::pair<std::string, std::string>> StrPairVec;

typedef std::shared_ptr<std::string> StrSPtr;
typedef std::shared_ptr<std::vector<std::string>> StrVecSPtr;

typedef std::array<std::string, 3> StrArray3;
typedef std::shared_ptr<StrArray3> StrArray3SPtr;
typedef std::vector<StrArray3SPtr> StrArray3SPtrVec;

typedef struct kstring_t KString;
