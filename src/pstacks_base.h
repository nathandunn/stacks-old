#include <vector>
#include <map>

#include "DNANSeq.h"
#include "stacks.h"
#include "mstack.h"

int    call_consensus(std::map<int, MergedStack *>& merged, std::map<int, PStack *>& unique, bool invoke_model);
int    call_alleles(MergedStack* mstack, std::vector<DNANSeq *>& reads);
int    write_results(map<int, MergedStack *> &, map<int, PStack *> &);
