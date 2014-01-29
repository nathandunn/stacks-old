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
@class ColorGenerator;

@interface DatumMO : NSManagedObject

@property (nonatomic, copy) NSString * name;
@property (nonatomic, copy) NSNumber * sampleId;
@property (nonatomic, copy) NSNumber * tagId;
//@property (nonatomic, copy) NSSet *alleles;
@property (nonatomic, copy) NSString * alleleData;
//@property (nonatomic, copy) ConsensusStackEntryMO *consensus;
//@property (nonatomic, copy) NSSet *depths;
@property (nonatomic, copy) NSString * depthData;
//@property (nonatomic, copy) NSSet *haplotypes;
@property (nonatomic, copy) NSData * haplotypeData;
@property (nonatomic, retain) LocusMO *locus;
//@property (nonatomic, retain) ModelStackEntryMO *model;
//@property (nonatomic, retain) ReferenceStackEntryMO *reference;
@property (nonatomic, retain) SampleMO *sample;
//@property (nonatomic, retain) NSSet *snps;
//@property (nonatomic, copy) NSString *snpData;
@property (nonatomic, copy) NSData *snpData;
//@property (nonatomic, copy) NSSet *stackEntries;
@property (nonatomic, copy) NSData *stackData;
//@property (nonatomic, copy) NSData * haploytpeData;

@property (nonatomic, retain) ColorGenerator *colorGenerator;
@end

@interface DatumMO (CoreDataGeneratedAccessors)


- (NSDictionary*)generateColorForOrder:(NSUInteger) order;

//- (void)addAllelesObject:(DatumAlleleMO *)value;
//- (void)removeAllelesObject:(DatumAlleleMO *)value;
//- (void)addAlleles:(NSSet *)values;
//- (void)removeAlleles:(NSSet *)values;
//
//- (void)addDepthsObject:(DepthMO *)value;
//- (void)removeDepthsObject:(DepthMO *)value;
//- (void)addDepths:(NSSet *)values;
//- (void)removeDepths:(NSSet *)values;
//
//- (void)addHaplotypesObject:(HaplotypeMO *)value;
//- (void)removeHaplotypesObject:(HaplotypeMO *)value;
//- (void)addHaplotypes:(NSSet *)values;
//- (void)removeHaplotypes:(NSSet *)values;
//
//- (void)addSnpsObject:(DatumSnpMO *)value;
//- (void)removeSnpsObject:(DatumSnpMO *)value;
//- (void)addSnps:(NSSet *)values;
//- (void)removeSnps:(NSSet *)values;

//- (void)addStackEntriesObject:(StackEntryMO *)value;
//- (void)removeStackEntriesObject:(StackEntryMO *)value;
//- (void)addStackEntries:(NSSet *)values;
//- (void)removeStackEntries:(NSSet *)values;


//- (NSArray*)getOrderedStackEntries;

@end
