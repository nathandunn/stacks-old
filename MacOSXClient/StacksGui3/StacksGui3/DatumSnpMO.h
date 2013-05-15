//
//  DatumSnpMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>
#import "SnpMO.h"

@class DatumMO;

@interface DatumSnpMO : SnpMO

@property (nonatomic, retain) DatumMO *datum;

@end
