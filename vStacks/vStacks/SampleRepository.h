//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//
//
//


#import <Foundation/Foundation.h>

@class SampleMO;


@interface SampleRepository : NSObject

+ (SampleRepository *)sharedInstance;
- (SampleMO *)getSampleForName:(NSString *)string andContext:(NSManagedObjectContext *)context andError:(NSError*)error;

- (SampleMO *)insertSample:(NSManagedObjectContext *)context id:(NSNumber *)number name:(NSString *)name;

- (NSArray *)getAllSamples:(NSManagedObjectContext *)context ;
@end