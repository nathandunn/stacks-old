//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StackRepository.h"
#import "StackMO.h"
#import "DatumMO.h"


@implementation StackRepository {

}
- (StackMO *)insertStack:(NSManagedObjectContext *)context datum:(DatumMO *)datum {
    StackMO *stackMO = [NSEntityDescription insertNewObjectForEntityForName:@"Stack" inManagedObjectContext:context];
    stackMO.datum = datum;
    datum.stack = stackMO;
    return stackMO;
}
@end