//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class SampleMO;


@interface SampleRepository : NSObject

+ (SampleRepository *)sharedInstance;
- (SampleMO *)getSampleForName:(NSString *)string andContext:(NSManagedObjectContext *)context andError:(NSError*)error;

//- (SampleMO *)insertSampleWithId:(NSNumber *)number andName:(NSString *)name;

- (SampleMO *)insertSample:(NSManagedObjectContext *)context id:(NSNumber *)number name:(NSString *)name;

- (NSArray *)getAllSamples:(NSManagedObjectContext *)context ;
@end