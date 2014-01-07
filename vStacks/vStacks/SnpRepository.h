//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class LocusSnpMO ;
@class LocusMO ;
@class DatumSnpMO ;
@class DatumMO ;

@interface SnpRepository : NSObject
+ (SnpRepository *)sharedInstance;
- (LocusSnpMO *)insertLocusSnp:(NSManagedObjectContext *)context column:(NSNumber *)column lratio:(NSNumber *)lratio rank1:(NSNumber *)rank1 rank2:(NSNumber *)rank2 rank3:(NSNumber *)rank3 rank4:(NSNumber *)rank4 locus:(LocusMO *)locus;
- (DatumSnpMO *)insertDatumSnp:(NSManagedObjectContext *)context column:(NSNumber *)column lratio:(NSNumber *)lratio rank1:(NSNumber *)rank1 rank2:(NSNumber *)rank2 rank3:(NSNumber *)rank3 rank4:(NSNumber *)rank4 datum:(DatumMO *)datum;
@end