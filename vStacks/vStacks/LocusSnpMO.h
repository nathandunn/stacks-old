//
//  LocusSnpMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>
#import "SnpMO.h"

@class LocusMO;

@interface LocusSnpMO : SnpMO

@property (nonatomic, retain) LocusMO *locus;

@end
