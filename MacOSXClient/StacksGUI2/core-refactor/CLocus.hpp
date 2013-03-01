//
//  CLocus.h
//  StacksGUI2
//
//  Created by Nathan Dunn on 3/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#ifndef StacksGUI2_CLocus_h
#define StacksGUI2_CLocus_h

#include <string>
using std::string;
#include <map>
using std::map;

#include "stacks.h"
//
// Catalog Locus Class
//
class CLocus : public Locus {
public:
    CLocus() : Locus() { this->f = 0.0; this->hcnt = 0; this->gcnt = 0; this->trans_gcnt = 0; };
    string annotation;
    string marker;
    double f;                 // Inbreeder's coefficient
    map<string, string> gmap; // Observed haplotype to genotype map for this locus.
    int hcnt;                 // Number of progeny containing a haplotype for this locus.
    int gcnt;                 // Number of progeny containing a valid genotype.
    int trans_gcnt;           // Number of progeny containing a valid 
                              // genotype, translated for a particular map type.
};

#endif
