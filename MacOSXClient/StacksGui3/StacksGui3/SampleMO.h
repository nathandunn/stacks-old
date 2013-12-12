//
//  SampleMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 12/12/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DatumMO, MatchMO, PopulationMO;

@interface SampleMO : NSManagedObject

@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) NSSet *datums;
@property (nonatomic, retain) PopulationMO *population;
@property (nonatomic, retain) NSSet *matches;
@end

@interface SampleMO (CoreDataGeneratedAccessors)

- (void)addDatumsObject:(DatumMO *)value;
- (void)removeDatumsObject:(DatumMO *)value;
- (void)addDatums:(NSSet *)values;
- (void)removeDatums:(NSSet *)values;

- (void)addMatchesObject:(MatchMO *)value;
- (void)removeMatchesObject:(MatchMO *)value;
- (void)addMatches:(NSSet *)values;
- (void)removeMatches:(NSSet *)values;

@end
