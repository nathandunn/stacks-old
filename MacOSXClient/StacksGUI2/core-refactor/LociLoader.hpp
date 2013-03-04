//
// Created by ndunn on 3/4/13.
//
// To change the template use AppCode | Preferences | File Templates.
//



#ifndef __LociLoader_H_
#define __LociLoader_H_

//double minor_allele_freq = 0;
//int progeny_limit = 0;

//#include <iostream>
#include "locus.h"
#include "PopMap.h"


/**
* This is code largely stolen from populations.cc
*/
class LociLoader {

public:
    int tabulate_haplotypes(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap);

    int create_genotype_map(CSLocus *pLocus, PopMap<CSLocus> *pMap);

    int call_population_genotypes(CSLocus *pLocus, PopMap<CSLocus> *pMap);



};

#endif //__LociLoader_H_
