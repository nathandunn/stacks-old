//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "LocusRepository.h"
#import "LocusMO.h"


@implementation LocusRepository {

}
- (LocusMO *)insertNewLocus:(NSManagedObjectContext *)model withId:(NSNumber *)id andConsensus:(NSString *)consensus andMarker:(NSString *)markers {
    LocusMO *locusMO = [NSEntityDescription insertNewObjectForEntityForName:@"Locus" inManagedObjectContext:model];
    locusMO.locusId = id ;
    locusMO.consensus = consensus;
    locusMO.marker = markers;
    return locusMO ;
}
@end