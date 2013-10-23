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
    
    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"chromosome!=nil"];
    [request1 setPredicate:predicate1];
    NSError *error1;
    NSArray *locusArray = [context executeFetchRequest:request1 error:&error1];
    return locusArray ;
}

@end