//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class PopulationMO;


@interface PopulationRepository : NSObject
- (PopulationMO *)insertPopulation:(NSManagedObjectContext *)context id:(NSNumber *)id name:(NSString *)name;

- (NSArray *)getAllPopulations:(NSManagedObjectContext *)context;

- (PopulationMO *)getPopulation:(NSManagedObjectContext *)context name:(NSString *)name;
@end