#pragma once
enum ReadsGroupBy { UMI = 1, POS = 2, READ_GROUP = 4, BARCODE = 8 };
enum DedupStrategy { RANDOM = 1, UNIQUE = 2, LONGEST = 4, MERGE = 8 };
