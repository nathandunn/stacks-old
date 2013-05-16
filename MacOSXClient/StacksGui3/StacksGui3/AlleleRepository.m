//
// Created by Nathan Dunn on 5/15/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "AlleleRepository.h"
#import "DatumAlleleMO.h"
#import "DatumMO.h"


@implementation AlleleRepository {

}
- (DatumAlleleMO *)insertDatumAllele:(NSManagedObjectContext *)context ratio:(NSNumber *)ratio depth:(NSNumber *)depth allele:(NSNumber *)allele datum:(DatumMO *)datum {
    DatumAlleleMO *datumAlleleMO = [NSEntityDescription insertNewObjectForEntityForName:@"DatumAllele" inManagedObjectContext:context];

    datumAlleleMO.depth = depth ;
    datumAlleleMO.ratio = ratio ;
    datumAlleleMO.allele = allele;
    datumAlleleMO.datum = datum ;

    return datumAlleleMO;
}
@end