//
//  DepthMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class GenotypeMO;

@interface DepthMO : NSManagedObject

@property (nonatomic, retain) NSNumber * depth;
@property (nonatomic, retain) GenotypeMO *genotype;

@end
