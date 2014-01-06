//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "LocusRepository.h"
#import "LocusMO.h"


@implementation LocusRepository {

}

+ (LocusRepository *)sharedInstance
{
    static LocusRepository *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[LocusRepository alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}

- (LocusMO *)insertNewLocus:(NSManagedObjectContext *)model withId:(NSNumber *)id andConsensus:(NSString *)consensus andMarker:(NSString *)markers {
    LocusMO *locusMO = [NSEntityDescription insertNewObjectForEntityForName:@"Locus" inManagedObjectContext:model];
    locusMO.locusId = id ;
    locusMO.consensus = consensus;
    locusMO.marker = markers;
    return locusMO ;
}

- (NSArray *)getAllLoci:(NSManagedObjectContext *)context {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Locus" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];
    NSError *error1;
    return [context executeFetchRequest:request1 error:&error1];
}

- (LocusMO *)getLocus:(NSManagedObjectContext *)context forId:(NSInteger)id {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Locus" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];

    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"locusId == %ld ", id];
    [request1 setPredicate:predicate1];
    NSError *error1;
    NSArray *locusArray = [context executeFetchRequest:request1 error:&error1];
    if(locusArray!=nil && locusArray.count==1){
        return [locusArray objectAtIndex:0] ;
    }
    else{
        return nil ;
    }
}

- (NSArray *)getLociWithChromsomes:(NSManagedObjectContext *)context{
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Locus" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];

    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"chromosome != nil"];
    NSPredicate *predicate2 = [NSPredicate predicateWithFormat:@"chromosome != ''"];
    [request1 setPredicate:predicate1];
    [request1 setPredicate:predicate2];
    NSError *error1;
    NSArray *locusArray = [context executeFetchRequest:request1 error:&error1];
    return locusArray ;
}

- (NSArray *)getLociLocations:(NSManagedObjectContext *)context{
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Locus" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];

    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"chromosome!=nil"];
    NSSortDescriptor *sortDescriptor = [[NSSortDescriptor alloc] initWithKey:@"chromosome" ascending:YES];
    NSArray *sortDescriptors = [[NSArray alloc] initWithObjects:sortDescriptor, nil];
    [request1 setPredicate:predicate1];
    [request1 setSortDescriptors:sortDescriptors];

    request1.resultType = NSDictionaryResultType;
    request1.propertiesToFetch = [NSArray arrayWithObject:[[entityDescription1 propertiesByName] objectForKey:@"chromosome"]];
    request1.returnsDistinctResults = YES;

    NSError *error1;
    NSArray *locusArray = [context executeFetchRequest:request1 error:&error1];


    NSMutableArray *returnArray = [[NSMutableArray alloc] init];
    for(NSDictionary *l in locusArray){
//        NSLog(@"out %@",l);
//        NSLog(@"out2 %@",l.allKeys);
        for(id key in l.allKeys){
//            NSLog(@"out3 %@", [l valueForKey:key]);
            [returnArray addObject:[l valueForKey:key]];
        }
    }

    return returnArray;
}

- (double)getMaxLocation:(NSManagedObjectContext *)context {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Locus" inManagedObjectContext:context];
    NSFetchRequest *request = [[NSFetchRequest alloc] init];
    NSSortDescriptor *sortDescriptor = [[NSSortDescriptor alloc] initWithKey:@"basePairs" ascending:NO];
    NSArray *sortDescriptors = [[NSArray alloc] initWithObjects:sortDescriptor, nil];


    [request setEntity:entityDescription1];
    [request setFetchLimit:1];
    [request setSortDescriptors:sortDescriptors];



    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"chromosome!=nil"];
    [request setPredicate:predicate1];
    NSError *error1;
    NSArray *locusArray = [context executeFetchRequest:request error:&error1];
    if(locusArray!=nil && locusArray.count==1){
        double doubleValue = [[[locusArray objectAtIndex:0] basePairs] doubleValue];
        NSLog(@"max value is %f",doubleValue);
        return doubleValue ;
    }
    return -1.0;
}
@end