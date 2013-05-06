//
//  SampleMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/6/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DatumMO, PopulationMO;

@interface SampleMO : NSManagedObject

@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) PopulationMO *population;
@property (nonatomic, retain) NSSet *datums;
@end

@interface SampleMO (CoreDataGeneratedAccessors)

- (void)addDatumsObject:(DatumMO *)value;
- (void)removeDatumsObject:(DatumMO *)value;
- (void)addDatums:(NSSet *)values;
- (void)removeDatums:(NSSet *)values;

@end
