//
//  Depth.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class Genotype;

@interface Depth : NSManagedObject

@property (nonatomic, retain) NSNumber * depth;
@property (nonatomic, retain) Genotype *genotype;

@end
