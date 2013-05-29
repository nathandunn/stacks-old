//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "PopulationRepository.h"
#import "PopulationMO.h"


@implementation PopulationRepository {

}

- (PopulationMO *)insertPopulation:(NSManagedObjectContext *)context id:(NSNumber *)id name:(NSString *)name {
    PopulationMO *populationMO = [NSEntityDescription insertNewObjectForEntityForName:@"Population" inManagedObjectContext:context];
    populationMO.populationId = id;
    populationMO.name = name ;
    return populationMO;
}

- (NSArray *)getAllPopulations:(NSManagedObjectContext *)context {
    NSEntityDescription *entityDescription2 = [NSEntityDescription
            entityForName:@"Population" inManagedObjectContext:context];
    NSFetchRequest *request2 = [[NSFetchRequest alloc] init];
    [request2 setEntity:entityDescription2];

    NSError *error;
    NSArray *populationArray = [context executeFetchRequest:request2 error:&error];
    return populationArray;
}

- (PopulationMO *)getPopulation:(NSManagedObjectContext *)context name:(NSString *)populationName{
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Population" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"name == %@", populationName];
    [request1 setPredicate:predicate1];
    [request1 setEntity:entityDescription1];
    NSError *error1;
    NSArray *populationArray = [context executeFetchRequest:request1 error:&error1];
    if (populationArray == nil || populationArray.count == 0) {
        PopulationMO *newPopulationMO = [NSEntityDescription insertNewObjectForEntityForName:@"Population" inManagedObjectContext:context];
        newPopulationMO.name = populationName;
        return newPopulationMO ;
    }
    else{
        return nil ;
    }
}

- (PopulationMO *)getPopulation:(NSManagedObjectContext *)context byIndexSortedByName:(NSInteger)index {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Population" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];
    NSError *error1;
    NSArray *populationArray = [context executeFetchRequest:request1 error:&error1];
    if (populationArray != nil || populationArray.count > 0) {
        NSSortDescriptor *descriptor = [NSSortDescriptor sortDescriptorWithKey:@"populationId"
                                                                     ascending:YES];
        NSArray *sortedArray = [populationArray sortedArrayUsingDescriptors:[NSArray arrayWithObject:descriptor]];
        PopulationMO *returnedPopulation = [sortedArray objectAtIndex:index] ;
        return returnedPopulation ;
    }
    else{
        return nil ;
    }
}
@end