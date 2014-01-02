//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class DatumMO;
@class SampleMO ;
@class LocusMO;
@class PopulationMO;

@interface DatumRepository : NSObject
+ (DatumRepository *)sharedInstance;

- (DatumMO *)insertDatum:(NSManagedObjectContext *)context name:(NSString *)name sampleId:(NSNumber *)id sample:(SampleMO *)sample locus:(LocusMO *)locus;

- (DatumMO *)getDatum:(NSManagedObjectContext *)context locusId:(NSInteger)locusId andSampleId:(NSInteger)sampleId;

- (NSArray *)getDatums:(NSManagedObjectContext *)context locus:(NSNumber *)locus andPopulation:(PopulationMO *)population;

- (NSArray *)getAllDatum:(NSManagedObjectContext *)context;

- (NSArray *)getDatumsOrdered:(NSManagedObjectContext *)context locus:(NSNumber *)locus andPopulation:(PopulationMO *)population;
@end