//
//  DatumMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class AlleleMO, DepthMO, HaplotypeMO, LocusMO, SampleMO, SnpMO, StackEntryMO;

@interface DatumMO : NSManagedObject

@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) NSNumber * tagId;
@property (nonatomic, retain) NSSet *alleles;
@property (nonatomic, retain) NSSet *depths;
@property (nonatomic, retain) NSSet *haplotypes;
@property (nonatomic, retain) LocusMO *locus;
@property (nonatomic, retain) SampleMO *sample;
@property (nonatomic, retain) NSSet *snps;
@property (nonatomic, retain) StackEntryMO *consensus;
@property (nonatomic, retain) StackEntryMO *model;
@property (nonatomic, retain) StackEntryMO *reference;
@property (nonatomic, retain) NSSet *stackEntries;
@end

@interface DatumMO (CoreDataGeneratedAccessors)

- (void)addAllelesObject:(AlleleMO *)value;
- (void)removeAllelesObject:(AlleleMO *)value;
- (void)addAlleles:(NSSet *)values;
- (void)removeAlleles:(NSSet *)values;

- (void)addDepthsObject:(DepthMO *)value;
- (void)removeDepthsObject:(DepthMO *)value;
- (void)addDepths:(NSSet *)values;
- (void)removeDepths:(NSSet *)values;

- (void)addHaplotypesObject:(HaplotypeMO *)value;
- (void)removeHaplotypesObject:(HaplotypeMO *)value;
- (void)addHaplotypes:(NSSet *)values;
- (void)removeHaplotypes:(NSSet *)values;

- (void)addSnpsObject:(SnpMO *)value;
- (void)removeSnpsObject:(SnpMO *)value;
- (void)addSnps:(NSSet *)values;
- (void)removeSnps:(NSSet *)values;

- (void)addStackEntriesObject:(StackEntryMO *)value;
- (void)removeStackEntriesObject:(StackEntryMO *)value;
- (void)addStackEntries:(NSSet *)values;
- (void)removeStackEntries:(NSSet *)values;

@end
