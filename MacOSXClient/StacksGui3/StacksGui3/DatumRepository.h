//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class DatumMO;
@class SampleMO ;

@interface DatumRepository : NSObject
- (DatumMO *)insertDatum:(NSManagedObjectContext *)context name:(NSString *)name sampleId:(NSNumber *)id sample:(SampleMO *)sample;

- (DatumMO *)getDatum:(NSManagedObjectContext *)context locusId:(NSInteger)locusId andSampleName:(NSString *)sampleName;
@end