//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "DatumRepository.h"
#import "DatumMO.h"
#import "SampleMO.h"
#import "LocusMO.h"
#import "PopulationMO.h"
#import "LocusRepository.h"


@implementation DatumRepository {

}

+ (DatumRepository *)sharedInstance
{
    static DatumRepository *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[DatumRepository alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}

- (DatumMO *)insertDatum:(NSManagedObjectContext *)context name:(NSString *)name sampleId:(NSNumber *)id sample:(SampleMO *)sample locus:(LocusMO *)locus {
    DatumMO *newDatumMO = [NSEntityDescription insertNewObjectForEntityForName:@"Datum" inManagedObjectContext:context];
    newDatumMO.name = name;
    newDatumMO.sampleId = id ;
    newDatumMO.sample = sample ;
    newDatumMO.locus = locus ;
    return newDatumMO ;
}


//- (DatumMO *)getDatum:(NSManagedObjectContext *)moc locusId:(NSInteger)locusId andSampleName:(NSString *)sampleName {
    - (DatumMO *)getDatum:(NSManagedObjectContext *)moc locusId:(NSInteger)locusId andSampleId:(NSInteger)sampleId {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Datum" inManagedObjectContext:moc];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];

//    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"locus.locusId == %ld and sample.name == %ld ", locusId, sampleName];
    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"tagId == %ld and sampleId == %ld ", locusId, sampleId];
    [request1 setPredicate:predicate1];
    NSError *error1;
    NSArray *datumArray = [moc executeFetchRequest:request1 error:&error1];
    if(datumArray!=nil && datumArray.count==1){
        return [datumArray objectAtIndex:0] ;
    }
    else{
        return nil ;
    }
}

- (NSArray *)getDatums:(NSManagedObjectContext *)context locus:(NSNumber *)locus andPopulation:(PopulationMO *)population {

    LocusMO *locusMO = [[LocusRepository sharedInstance] getLocus:context forId:locus.integerValue];
//    NSLog(@"locus %@",locusMO);


    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Datum" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
//    NSLog(@"locusID: %@ population: %@",locus,population);
    NSMutableArray *sampleIds = [[NSMutableArray alloc] init];

    for(SampleMO *sampleMO in population.samples){
        [sampleIds addObject:sampleMO.sampleId];
    }
    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"tagId == %@ and sampleId in (%@) ", locus,sampleIds];
    [request1 setPredicate:predicate1];
    [request1 setEntity:entityDescription1];
    NSError *error1;
    return [context executeFetchRequest:request1 error:&error1];
}

- (NSArray *)getAllDatum:(NSManagedObjectContext *)context {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Datum" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];
    NSError *error1;
    return [context executeFetchRequest:request1 error:&error1];
}

- (NSArray *)getDatumsOrdered:(NSManagedObjectContext *)context locus:(NSNumber *)locus andPopulation:(PopulationMO *)population {
    NSArray *unsortedArray = [self getDatums:context locus:locus andPopulation:population];
//    NSLog(@"datums returned: %ld",unsortedArray.count) ;
    NSArray *sortedArray;
    NSNumberFormatter *numberFormatter = [[NSNumberFormatter alloc] init];
    numberFormatter.numberStyle = NSNumberFormatterNoStyle;

    sortedArray = [unsortedArray sortedArrayUsingComparator:^NSComparisonResult(id a, id b) {
        NSString *first = [(DatumMO*)a name];
        NSString *second = [(DatumMO*)b name];

        int firstScore = 0 ;
        int secondScore = 0 ;



        if([first isEqualToString:@"male"]) firstScore -= 1000 ;
        if([second isEqualToString:@"male"]) secondScore -= 1000 ;
        if([first isEqualToString:@"female"]) firstScore -= 100 ;
        if([second isEqualToString:@"female"]) secondScore -= 100 ;

        // handle sample_ vs progeny_

        NSArray *stringOne  = [first componentsSeparatedByString:@"_"];
        NSArray *stringTwo = [second componentsSeparatedByString:@"_"];

        if(stringOne.count>1){
            firstScore += [[numberFormatter numberFromString:[stringOne objectAtIndex:1]] intValue];
        }
        if(stringTwo.count>1){
            secondScore += [[numberFormatter numberFromString:[stringTwo objectAtIndex:1]] intValue];
        }


//        if([first isEqualToString:@"female"] && [second isEqualToString:@"male"]){
//          return NSOrderedDescending;
//        }
//        if([first isEqualToString:@"male"] && [second isEqualToString:@"female"]){
//            return NSOrderedAscending;
//        }

        return (firstScore-secondScore<0) ? NSOrderedAscending : NSOrderedDescending ;
//
//        return [first compare:second];
    }];
    return sortedArray;
}

- (NSDictionary *)getDatums:(NSManagedObjectContext *)context forSample:(NSNumber*)sampleId{
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Datum" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];

    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@" sampleId == %@ ", sampleId];
    [request1 setPredicate:predicate1];
    [request1 setEntity:entityDescription1];
    NSError *error1;
    NSArray *datums = [context executeFetchRequest:request1 error:&error1];

//    NSMutableDictionary *dictionary = [[NSMutableDictionary alloc] initWithCapacity:datums.count];
    NSMutableDictionary *dictionary = [NSMutableDictionary dictionaryWithCapacity:datums.count];

    for(DatumMO *datumMO in datums){
//        [dictionary insertValue:datumMO inPropertyWithKey:datumMO.tagId.stringValue];
        [dictionary setObject:datumMO forKey:datumMO.tagId.stringValue];
    }

    return dictionary ;

}
@end