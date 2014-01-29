//
// Created by Nathan Dunn on 5/15/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "AlleleRepository.h"
#import "DatumAlleleMO.h"
#import "DatumMO.h"
#import "LocusAlleleMO.h"
#import "LocusMO.h"


@implementation AlleleRepository {

}

+ (AlleleRepository *)sharedInstance
{
    static AlleleRepository *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[AlleleRepository alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}

- (DatumAlleleMO *)insertDatumAllele:(NSManagedObjectContext *)context ratio:(NSNumber *)ratio depth:(NSNumber *)depth allele:(NSNumber *)allele datum:(DatumMO *)datum {
    DatumAlleleMO *datumAlleleMO = [NSEntityDescription insertNewObjectForEntityForName:@"DatumAllele" inManagedObjectContext:context];

    datumAlleleMO.depth = depth ;
    datumAlleleMO.ratio = ratio ;
    datumAlleleMO.allele = allele;
    datumAlleleMO.datum = datum ;

    return datumAlleleMO;
}

- (LocusAlleleMO *)insertLocusAllele:(NSManagedObjectContext *)context depth:(NSNumber *)depth allele:(NSNumber *)allele locus:(LocusMO *)locus {
    LocusAlleleMO *locusAlleleMO = [NSEntityDescription insertNewObjectForEntityForName:@"LocusAllele" inManagedObjectContext:context];

    locusAlleleMO.depth = depth ;
    locusAlleleMO.allele = allele;
    locusAlleleMO.locus = locus ;

    return locusAlleleMO;
}
@end