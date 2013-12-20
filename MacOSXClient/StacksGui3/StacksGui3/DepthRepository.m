//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "DepthRepository.h"
#import "DepthMO.h"


@implementation DepthRepository {

}
+ (DepthRepository *)sharedInstance
{
    static DepthRepository *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[DepthRepository alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}


- (DepthMO *)insertDepth:(NSManagedObjectContext *)context depth:(NSNumber *)depth andOrder:(int)order {
    DepthMO *depthMO = [NSEntityDescription insertNewObjectForEntityForName:@"Depth" inManagedObjectContext:context];
    depthMO.depth = depth ;
    depthMO.order = [NSNumber numberWithInt:order];
    return depthMO ;
}
@end