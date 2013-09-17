//
//  HaplotypeMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 9/17/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DatumMO;

@interface HaplotypeMO : NSManagedObject

@property (nonatomic, retain) NSString * haplotype;
@property (nonatomic, retain) NSNumber * order;
@property (nonatomic, retain) DatumMO *datum;

@end
