//
//  LocusMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/6/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class AlleleMO, DatumMO, SampleMO, SnpMO;

@interface LocusMO : NSManagedObject

@property (nonatomic, retain) NSString * consensus;
@property (nonatomic, retain) NSNumber * length;
@property (nonatomic, retain) NSNumber * locusId;
@property (nonatomic, retain) NSString * marker;
@property (nonatomic, retain) NSString * ratio;
@property (nonatomic, retain) NSSet *alleles;
@property (nonatomic, retain) NSSet *datums;
@property (nonatomic, retain) NSSet *snps;
@property (nonatomic, retain) SampleMO *sample;
@end

@interface LocusMO (CoreDataGeneratedAccessors)

- (void)addAllelesObject:(AlleleMO *)value;
- (void)removeAllelesObject:(AlleleMO *)value;
- (void)addAlleles:(NSSet *)values;
- (void)removeAlleles:(NSSet *)values;

- (void)addDatumsObject:(DatumMO *)value;
- (void)removeDatumsObject:(DatumMO *)value;
- (void)addDatums:(NSSet *)values;
- (void)removeDatums:(NSSet *)values;

- (void)addSnpsObject:(SnpMO *)value;
- (void)removeSnpsObject:(SnpMO *)value;
- (void)addSnps:(NSSet *)values;
- (void)removeSnps:(NSSet *)values;

@end
