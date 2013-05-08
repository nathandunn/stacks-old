//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StackEntryRepository.h"
#import "StackEntryMO.h"
#import "StackMO.h"


@implementation StackEntryRepository {

}
- (StackEntryMO *)insertStackEntry:(NSManagedObjectContext *)context entryId:(NSNumber *)id relationship:(NSString *)relationship block:(NSString *)block sequenceId:(NSString *)sequenceId sequence:(NSString *)sequence stack:(StackMO *)stack {
    StackEntryMO *stackEntryMO = [NSEntityDescription insertNewObjectForEntityForName:@"StackEntry" inManagedObjectContext:context];
    stackEntryMO.entryId = id ;
    stackEntryMO.relationship = relationship;
    stackEntryMO.block = block ;
    stackEntryMO.sequenceId = sequenceId ;
    stackEntryMO.sequence = sequence ;
    stackEntryMO.stack = stack ;
    return stackEntryMO ;
}
@end