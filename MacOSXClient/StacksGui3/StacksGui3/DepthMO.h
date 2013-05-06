//
//  DepthMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/6/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DatumMO;

@interface DepthMO : NSManagedObject

@property (nonatomic, retain) NSNumber * depth;
@property (nonatomic, retain) DatumMO *datum;

@end
