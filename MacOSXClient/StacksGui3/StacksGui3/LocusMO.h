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

@property (nonatomic, retain) NSNumber * basePairs;
@property (nonatomic, retain) NSString * chromosome;
@property (nonatomic, retain) NSString * strand;
@property (nonatomic, retain) NSString * type;
@property (nonatomic, retain) NSString * consensus;
@property (nonatomic, retain) NSNumber * length;
@property (nonatomic, retain) NSNumber * locusId;
@property (nonatomic, retain) NSNumber * parentCount;
@property (nonatomic, retain) NSString * marker;
@property (nonatomic, retain) NSString * ratio;

//@property (nonatomic, retain) NSSet *alleles;
@property (nonatomic, retain) NSSet *datums;
//@property (nonatomic, retain) NSSet *snps;

@property (nonatomic, retain) NSString *alleleData;
//@property (nonatomic, retain) NSString *snpData;
@property (nonatomic, retain) NSData *snpData;

// not a managed part, populated by datum rendering . . .
@property (atomic, retain) NSDictionary *haplotypeOrder;
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
