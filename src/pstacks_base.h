#include <vector>
#include <map>

#include "DNANSeq.h"
#include "stacks.h"
#include "mstack.h"

extern string prefix_path;
extern int    sql_id;
extern modelt model_type;

int call_consensus(std::map<int, MergedStack *>& merged,
                   std::map<int, PStack *>& unique,
                   bool invoke_model);

int call_alleles(MergedStack* mstack,
                 std::vector<DNANSeq *>& reads);

int write_results(map<int, MergedStack*>& merged,
                  map<int, PStack*>& unique,
                  bool gzip, bool paired_end);
