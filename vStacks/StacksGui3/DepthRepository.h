//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class DepthMO;


@interface DepthRepository : NSObject

+ (DepthRepository *)sharedInstance;
- (DepthMO *)insertDepth:(NSManagedObjectContext *)context depth:(NSNumber *)depth andOrder:(int)order;
@end