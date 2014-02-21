//
//  PopulationMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class SampleMO;

@interface PopulationMO : NSManagedObject

@property (nonatomic, copy) NSString * name;
@property (nonatomic, copy) NSNumber * populationId;
@property (nonatomic, copy) NSSet *samples;
@property (nonatomic, copy) NSData *metaData;

- (NSString *)annotatedName;
@end

@interface PopulationMO (CoreDataGeneratedAccessors)

- (void)addSamplesObject:(SampleMO *)value;
- (void)removeSamplesObject:(SampleMO *)value;
- (void)addSamples:(NSSet *)values;
- (void)removeSamples:(NSSet *)values;

@end
