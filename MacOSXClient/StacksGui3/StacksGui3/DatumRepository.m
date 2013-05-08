//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "DatumRepository.h"
#import "DatumMO.h"
#import "SampleMO.h"


@implementation DatumRepository {

}
- (DatumMO *)insertDatum:(NSManagedObjectContext *)context name:(NSString *)name sampleId:(NSNumber *)id sample:(SampleMO *)sample {
    DatumMO *newDatumMO = [NSEntityDescription insertNewObjectForEntityForName:@"Datum" inManagedObjectContext:context];
    newDatumMO.name = name;
    newDatumMO.sampleId = id ;
    newDatumMO.sample = sample ;

    return newDatumMO ;
}


- (DatumMO *)getDatum:(NSManagedObjectContext *)moc locusId:(NSInteger)locusId andSampleName:(NSString *)sampleName {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Datum" inManagedObjectContext:moc];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"locus.locusId == %ld and sample.name == %@ ", locusId, sampleName];
    [request1 setPredicate:predicate1];
    [request1 setEntity:entityDescription1];
    NSError *error1;
    NSArray *datumArray = [moc executeFetchRequest:request1 error:&error1];
    if(datumArray!=nil && datumArray.count==1){
        return [datumArray objectAtIndex:0] ;
    }
    else{
        return nil ;
    }
}
@end