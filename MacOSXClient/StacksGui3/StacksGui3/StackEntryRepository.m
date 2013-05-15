//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StackEntryRepository.h"
#import "StackEntryMO.h"
#import "ConsensusStackEntryMO.h"
#import "DatumMO.h"
#import "ModelStackEntryMO.h"


@implementation StackEntryRepository {

}
- (StackEntryMO *)insertStackEntry:(NSManagedObjectContext *)context entryId:(NSNumber *)entryId relationship:(NSString *)relationship block:(NSString *)block sequenceId:(NSString *)sequenceId sequence:(NSString *)sequence datum:(DatumMO *)datum {
    StackEntryMO *stackEntryMO = [NSEntityDescription insertNewObjectForEntityForName:@"StackEntry" inManagedObjectContext:context];
    stackEntryMO.entryId = entryId;
    stackEntryMO.relationship = relationship;
    stackEntryMO.block = block ;
    stackEntryMO.sequenceId = sequenceId ;
    stackEntryMO.sequence = sequence ;
    stackEntryMO.datum = datum ;
    return stackEntryMO ;
}

- (ConsensusStackEntryMO *)insertConsensusStackEntry:(NSManagedObjectContext *)context entryId:(NSNumber *)entryId block:(NSString*)block sequenceId:(NSString*)sequenceId sequence:(NSString*)sequence datum:(DatumMO *)datum {
    ConsensusStackEntryMO *stackEntryMO = [NSEntityDescription insertNewObjectForEntityForName:@"ConsensusStackEntry" inManagedObjectContext:context];
    stackEntryMO.entryId = entryId;
    stackEntryMO.relationship = @"consensus" ;
    stackEntryMO.block = block ;
    stackEntryMO.sequenceId = sequenceId ;
    stackEntryMO.sequence = sequence ;
    stackEntryMO.datum = datum ;
    return stackEntryMO ;
}

- (ModelStackEntryMO *)insertModelStackEntry:(NSManagedObjectContext *)context entryId:(NSNumber *)entryId block:(NSString *)block sequenceId:(NSString *)sequenceId sequence:(id)sequence datum:(DatumMO *)datum {
    ModelStackEntryMO *stackEntryMO = [NSEntityDescription insertNewObjectForEntityForName:@"ModelStackEntry" inManagedObjectContext:context];
    stackEntryMO.entryId = entryId;
    stackEntryMO.relationship = @"model" ;
    stackEntryMO.block = block ;
    stackEntryMO.sequenceId = sequenceId ;
    stackEntryMO.sequence = sequence ;
    stackEntryMO.datum = datum ;
    return stackEntryMO ;
}
@end