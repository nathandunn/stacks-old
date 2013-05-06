//
//  StackMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/6/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DatumMO, StackEntryMO;

@interface StackMO : NSManagedObject

@property (nonatomic, retain) StackEntryMO *consensus;
@property (nonatomic, retain) DatumMO *datum;
@property (nonatomic, retain) StackEntryMO *model;
@property (nonatomic, retain) StackEntryMO *reference;
@property (nonatomic, retain) NSSet *stackEntries;
@end

@interface StackMO (CoreDataGeneratedAccessors)

- (void)addStackEntriesObject:(StackEntryMO *)value;
- (void)removeStackEntriesObject:(StackEntryMO *)value;
- (void)addStackEntries:(NSSet *)values;
- (void)removeStackEntries:(NSSet *)values;

@end
