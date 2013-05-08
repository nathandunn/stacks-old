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

- (PopulationMO *)insertPopulation:(NSManagedObjectContext *)context withId:(NSNumber *)id andName:(NSString *)name {
    PopulationMO *populationMO = [NSEntityDescription insertNewObjectForEntityForName:@"Population" inManagedObjectContext:context];
    populationMO.populationId = id;
    populationMO.name = name ;
    return populationMO;
}
@end