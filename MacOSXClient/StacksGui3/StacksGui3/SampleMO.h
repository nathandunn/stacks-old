//
//  SampleMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DatumMO, PopulationMO;

@interface SampleMO : NSManagedObject

@property (nonatomic, copy) NSString * name;
@property (nonatomic, copy) NSNumber * sampleId;
@property (nonatomic, copy) NSSet *datums;
@property (nonatomic, retain) PopulationMO *population;
@end

@interface SampleMO (CoreDataGeneratedAccessors)

- (void)addDatumsObject:(DatumMO *)value;
- (void)removeDatumsObject:(DatumMO *)value;
- (void)addDatums:(NSSet *)values;
- (void)removeDatums:(NSSet *)values;

@end
