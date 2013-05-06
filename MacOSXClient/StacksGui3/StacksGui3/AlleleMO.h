//
//  AlleleMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/6/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class LocusMO;

@interface AlleleMO : NSManagedObject

@property (nonatomic, retain) NSString * type;
@property (nonatomic, retain) NSString * value;
@property (nonatomic, retain) LocusMO *locus;

@end
