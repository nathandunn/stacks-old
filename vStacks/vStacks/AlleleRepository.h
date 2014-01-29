//
// Created by Nathan Dunn on 5/15/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class DatumAlleleMO;
@class DatumMO;
@class LocusAlleleMO;
@class LocusMO;


@interface AlleleRepository : NSObject
+ (AlleleRepository *)sharedInstance;

- (DatumAlleleMO *)insertDatumAllele:(NSManagedObjectContext *)context ratio:(NSNumber *)ratio depth:(NSNumber *)depth allele:(NSNumber *)allele datum:(DatumMO *)datum;

- (LocusAlleleMO *)insertLocusAllele:(NSManagedObjectContext *)context depth:(NSNumber *)column allele:(NSNumber *)allele locus:(LocusMO *)locus;
@end