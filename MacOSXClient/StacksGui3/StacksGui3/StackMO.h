//
//  StackMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class GenotypeMO;

@interface StackMO : NSManagedObject

@property (nonatomic, retain) GenotypeMO *genotype;
@property (nonatomic, retain) NSSet *stackEntries;
@property (nonatomic, retain) NSManagedObject *consensus;
@property (nonatomic, retain) NSManagedObject *model;
@property (nonatomic, retain) NSManagedObject *reference;
@end

@interface StackMO (CoreDataGeneratedAccessors)

- (void)addStackEntriesObject:(NSManagedObject *)value;
- (void)removeStackEntriesObject:(NSManagedObject *)value;
- (void)addStackEntries:(NSSet *)values;
- (void)removeStackEntries:(NSSet *)values;

@end
