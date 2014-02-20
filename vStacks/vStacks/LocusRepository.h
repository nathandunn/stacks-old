//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//
//
//


#import <Foundation/Foundation.h>

@class LocusMO ;

@interface LocusRepository : NSObject

+ (LocusRepository *)sharedInstance;
- (LocusMO *)insertNewLocus:(NSManagedObjectContext *)model withId:(NSNumber *)id andConsensus:(NSString *)consensus andMarker:(NSString *)markers;

- (NSArray *)getAllLoci:(NSManagedObjectContext *)context;

- (LocusMO *)getLocus:(NSManagedObjectContext *)context forId:(NSInteger)id;
- (NSArray *)getLociWithChromsomes:(NSManagedObjectContext *)context;
- (NSArray *)getLociLocations:(NSManagedObjectContext *)context;

- (double)getMaxLocation:(NSManagedObjectContext *)context;

//- (NSNumber *)getProgenyCount:(NSManagedObjectContext *)context locus:(LocusMO *)locus;

- (NSDictionary *)getAggregateProgenyCount:(NSManagedObjectContext *)context;

@end