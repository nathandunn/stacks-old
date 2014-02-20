//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//
//
//


#import <Foundation/Foundation.h>

@class PopulationMO;


@interface PopulationRepository : NSObject

+ (PopulationRepository *)sharedInstance;
- (PopulationMO *)insertPopulation:(NSManagedObjectContext *)context id:(NSNumber *)id name:(NSString *)name;

- (NSArray *)getAllPopulations:(NSManagedObjectContext *)context;

- (PopulationMO *)getPopulationOrCreate:(NSManagedObjectContext *)context name:(NSString *)name;

- (PopulationMO *)getPopulation:(NSManagedObjectContext *)context byIndexSortedByName:(NSInteger)index;

- (PopulationMO *)getPopulation:(NSManagedObjectContext *)context name:(NSString *)name;
@end