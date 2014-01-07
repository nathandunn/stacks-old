//
//  DatumAlleleMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>
#import "AlleleMO.h"

@class DatumMO;

@interface DatumAlleleMO : AlleleMO

@property (nonatomic, retain) DatumMO *datum;

@end
