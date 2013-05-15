//
//  DatumMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class ConsensusStackEntryMO, DatumAlleleMO, DatumSnpMO, DepthMO, HaplotypeMO, LocusMO, ModelStackEntryMO, ReferenceStackEntryMO, SampleMO, StackEntryMO;

@interface DatumMO : NSManagedObject

@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) NSNumber * tagId;
@property (nonatomic, retain) NSSet *alleles;
@property (nonatomic, retain) ConsensusStackEntryMO *consensus;
@property (nonatomic, retain) NSSet *depths;
@property (nonatomic, retain) NSSet *haplotypes;
@property (nonatomic, retain) LocusMO *locus;
@property (nonatomic, retain) ModelStackEntryMO *model;
@property (nonatomic, retain) ReferenceStackEntryMO *reference;
@property (nonatomic, retain) SampleMO *sample;
@property (nonatomic, retain) NSSet *snps;
@property (nonatomic, retain) NSSet *stackEntries;
@end

@interface DatumMO (CoreDataGeneratedAccessors)

- (void)addAllelesObject:(DatumAlleleMO *)value;
- (void)removeAllelesObject:(DatumAlleleMO *)value;
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

- (void)addSnpsObject:(DatumSnpMO *)value;
- (void)removeSnpsObject:(DatumSnpMO *)value;
- (void)addSnps:(NSSet *)values;
- (void)removeSnps:(NSSet *)values;

- (void)addStackEntriesObject:(StackEntryMO *)value;
- (void)removeStackEntriesObject:(StackEntryMO *)value;
- (void)addStackEntries:(NSSet *)values;
- (void)removeStackEntries:(NSSet *)values;

@end
