//
//  DatumMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "DatumMO.h"
#import "ConsensusStackEntryMO.h"
#import "DatumAlleleMO.h"
#import "DatumSnpMO.h"
#import "DepthMO.h"
#import "HaplotypeMO.h"
#import "LocusMO.h"
#import "ModelStackEntryMO.h"
#import "ReferenceStackEntryMO.h"
#import "SampleMO.h"
#import "StackEntryMO.h"


@implementation DatumMO

@dynamic name;
@dynamic sampleId;
@dynamic tagId;
@dynamic alleles;
@dynamic consensus;
@dynamic depths;
@dynamic haplotypes;
@dynamic locus;
@dynamic model;
@dynamic reference;
@dynamic sample;
@dynamic snps;
@dynamic stackEntries;

@end
