//
//  PopulationMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/6/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class SampleMO;

@interface PopulationMO : NSManagedObject

@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) NSNumber * populationId;
@property (nonatomic, retain) NSSet *samples;
@end

@interface PopulationMO (CoreDataGeneratedAccessors)

- (void)addSamplesObject:(SampleMO *)value;
- (void)removeSamplesObject:(SampleMO *)value;
- (void)addSamples:(NSSet *)values;
- (void)removeSamples:(NSSet *)values;

@end
