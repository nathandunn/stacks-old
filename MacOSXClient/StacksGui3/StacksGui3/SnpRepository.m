//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "SnpRepository.h"
#import "SnpMO.h"
#import "LocusMO.h"
#import "DatumMO.h"
#import "LocusSnpMO.h"
#import "DatumSnpMO.h"
#import "DatumSnpMO.h"


@implementation SnpRepository {

}

+ (SnpRepository *)sharedInstance
{
    static SnpRepository *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[SnpRepository alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}

- (LocusSnpMO *)insertLocusSnp:(NSManagedObjectContext *)context column:(NSNumber *)column lratio:(NSNumber *)lratio rank1:(NSNumber *)rank1 rank2:(NSNumber *)rank2 rank3:(NSNumber *)rank3 rank4:(NSNumber *)rank4 locus:(LocusMO *)locus {
    LocusSnpMO *snpMO = [NSEntityDescription insertNewObjectForEntityForName:@"LocusSnp" inManagedObjectContext:context];

    snpMO.column = column ;
    snpMO.lratio = lratio ;
    snpMO.rank1 = rank1 ;
    snpMO.rank2 = rank2 ;
    snpMO.rank3 = rank3 ;
    snpMO.rank4 = rank4 ;
    snpMO.locus = locus ;

    return snpMO;
}


- (DatumSnpMO *)insertDatumSnp:(NSManagedObjectContext *)context column:(NSNumber *)column lratio:(NSNumber *)lratio rank1:(NSNumber *)rank1 rank2:(NSNumber *)rank2 rank3:(NSNumber *)rank3 rank4:(NSNumber *)rank4 datum:(DatumMO *)datum{
    DatumSnpMO *snpMO = [NSEntityDescription insertNewObjectForEntityForName:@"DatumSnp" inManagedObjectContext:context];

    snpMO.column = column ;
    snpMO.lratio = lratio ;
    snpMO.rank1 = rank1 ;
    snpMO.rank2 = rank2 ;
    if(rank3!=nil){
        snpMO.rank3 = rank3 ;
        snpMO.rank4 = rank4 ;
    }
    snpMO.datum = datum ;

    return snpMO;
}


@end