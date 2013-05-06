//
//  DatumMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/6/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class AlleleMO, DepthMO, HaplotypeMO, LocusMO, SnpMO, StackMO;

@interface DatumMO : NSManagedObject

@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) NSNumber * tagId;
@property (nonatomic, retain) NSSet *depths;
@property (nonatomic, retain) NSSet *haplotypes;
@property (nonatomic, retain) LocusMO *locus;
@property (nonatomic, retain) NSSet *stacks;
@property (nonatomic, retain) NSSet *alleles;
@property (nonatomic, retain) NSSet *snps;
@end

@interface DatumMO (CoreDataGeneratedAccessors)

- (void)addDepthsObject:(DepthMO *)value;
- (void)removeDepthsObject:(DepthMO *)value;
- (void)addDepths:(NSSet *)values;
- (void)removeDepths:(NSSet *)values;

- (void)addHaplotypesObject:(HaplotypeMO *)value;
- (void)removeHaplotypesObject:(HaplotypeMO *)value;
- (void)addHaplotypes:(NSSet *)values;
- (void)removeHaplotypes:(NSSet *)values;

- (void)addStacksObject:(StackMO *)value;
- (void)removeStacksObject:(StackMO *)value;
- (void)addStacks:(NSSet *)values;
- (void)removeStacks:(NSSet *)values;

- (void)addAllelesObject:(AlleleMO *)value;
- (void)removeAllelesObject:(AlleleMO *)value;
- (void)addAlleles:(NSSet *)values;
- (void)removeAlleles:(NSSet *)values;

- (void)addSnpsObject:(SnpMO *)value;
- (void)removeSnpsObject:(SnpMO *)value;
- (void)addSnps:(NSSet *)values;
- (void)removeSnps:(NSSet *)values;

@end
