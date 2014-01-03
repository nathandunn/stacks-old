//
//  LocusMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DatumMO, LocusAlleleMO, LocusSnpMO, DepthMO;
@class ColorGenerator;

@interface LocusMO : NSManagedObject

@property (nonatomic, copy) NSNumber * basePairs;
@property (nonatomic, copy) NSString * chromosome;
@property (nonatomic, copy) NSString * strand;
@property (nonatomic, copy) NSString * type;
@property (nonatomic, copy) NSString * consensus;
@property (nonatomic, copy) NSNumber * length;
@property (nonatomic, copy) NSNumber * locusId;
@property (nonatomic, copy) NSNumber * parentCount;
@property (nonatomic, copy) NSString * marker;
@property (nonatomic, copy) NSString * ratio;

//@property (nonatomic, copy) NSSet *alleles;
@property (nonatomic, strong) NSSet *datums;
//@property (nonatomic, copy) NSSet *snps;

@property (nonatomic, copy) NSString *alleleData;
//@property (nonatomic, copy) NSString *snpData;
@property (nonatomic, strong) NSData *snpData;

// not a managed part, populated by datum rendering . . .
@property (atomic, strong) NSDictionary *haplotypeOrder;
@property (atomic, retain) ColorGenerator *colorGenerator;


- (NSUInteger)lookupHaplotypeOrder:(NSString *)haplotype;

- (NSUInteger)lookupDepthOrder:(DepthMO *)mo;
@end

@interface LocusMO (CoreDataGeneratedAccessors)

//- (void)addAllelesObject:(LocusAlleleMO *)value;
//- (void)removeAllelesObject:(LocusAlleleMO *)value;
//- (void)addAlleles:(NSSet *)values;
//- (void)removeAlleles:(NSSet *)values;

- (void)addDatumsObject:(DatumMO *)value;
- (void)removeDatumsObject:(DatumMO *)value;
- (void)addDatums:(NSSet *)values;
- (void)removeDatums:(NSSet *)values;

//- (void)addSnpsObject:(LocusSnpMO *)value;
//- (void)removeSnpsObject:(LocusSnpMO *)value;
//- (void)addSnps:(NSSet *)values;
//- (void)removeSnps:(NSSet *)values;

@end
