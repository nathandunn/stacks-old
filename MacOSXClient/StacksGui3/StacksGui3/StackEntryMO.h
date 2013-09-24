//
//  StackEntryMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class DatumMO;

@interface StackEntryMO : NSManagedObject

@property (nonatomic, retain) NSString * block;
@property (nonatomic, retain) NSNumber * entryId;
@property (nonatomic, retain) NSString * relationship;
@property (nonatomic, retain) NSString * sequence;
@property (nonatomic, retain) NSString * sequenceId;
@property (nonatomic, retain) DatumMO *datum;

@end
