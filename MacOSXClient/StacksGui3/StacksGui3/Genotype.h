//
//  Genotype.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class Depth, Locus, Stack;

@interface Genotype : NSManagedObject

@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) NSNumber * tagId;
@property (nonatomic, retain) Locus *locus;
@property (nonatomic, retain) NSSet *stacks;
@property (nonatomic, retain) NSManagedObject *haplotypes;
@property (nonatomic, retain) Depth *depths;
@end

@interface Genotype (CoreDataGeneratedAccessors)

- (void)addStacksObject:(Stack *)value;
- (void)removeStacksObject:(Stack *)value;
- (void)addStacks:(NSSet *)values;
- (void)removeStacks:(NSSet *)values;

@end
