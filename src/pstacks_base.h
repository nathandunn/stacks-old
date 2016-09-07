#include <vector>
#include <map>

#include "DNANSeq.h"
#include "stacks.h"
#include "mstack.h"

int    populate_merged_tags(std::map<int, PStack *>& unique, std::map<int, MergedStack *>& merged);
int    call_consensus(std::map<int, MergedStack *>& merged, std::map<int, PStack *>& unique, bool invoke_model);
int    call_alleles(MergedStack* mstack, std::vector<DNANSeq *>& reads);
