#include <vector>
#include <map>

#include "models.h"
#include "DNANSeq.h"
#include "stacks.h"
#include "mstack.h"

extern string prefix_path;
extern int    sql_id;
extern modelt model_type;

int call_consensus(map<int, MergedStack *>& merged,
                   map<int, PStack *>& unique,
                   bool invoke_model);

void call_alleles(MergedStack* mstack,
                 vector<DNANSeq *>& reads);

void calc_coverage_distribution(const map<int, PStack*>& unique,
                                const map<int, MergedStack *>& merged,
                                double& mean,
                                double& stdev,
                                double& max);

int write_results(map<int, MergedStack*>& merged,
                  map<int, PStack*>& unique,
                  bool gzip, bool paired_end);
