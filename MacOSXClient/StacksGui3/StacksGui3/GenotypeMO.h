//
//  GenotypeMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DepthMO, LocusMO, StackMO;

@interface GenotypeMO : NSManagedObject

@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) NSNumber * tagId;
@property (nonatomic, retain) LocusMO *locus;
@property (nonatomic, retain) NSSet *stacks;
@property (nonatomic, retain) NSManagedObject *haplotypes;
@property (nonatomic, retain) DepthMO *depths;
@end

@interface GenotypeMO (CoreDataGeneratedAccessors)

- (void)addStacksObject:(StackMO *)value;
- (void)removeStacksObject:(StackMO *)value;
- (void)addStacks:(NSSet *)values;
- (void)removeStacks:(NSSet *)values;

@end
