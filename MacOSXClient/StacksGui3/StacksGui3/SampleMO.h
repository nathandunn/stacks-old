//
//  SampleMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/6/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class LocusMO, PopulationMO;

@interface SampleMO : NSManagedObject

@property (nonatomic, retain) PopulationMO *population;
@property (nonatomic, retain) NSSet *loci;
@end

@interface SampleMO (CoreDataGeneratedAccessors)

- (void)addLociObject:(LocusMO *)value;
- (void)removeLociObject:(LocusMO *)value;
- (void)addLoci:(NSSet *)values;
- (void)removeLoci:(NSSet *)values;

@end
