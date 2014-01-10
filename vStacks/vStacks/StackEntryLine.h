//
//  StackEntryLine.h
//  vStacks
//
//  Created by Nathan Dunn on 1/10/14.
//  Copyright (c) 2014 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>


@interface StackEntryLine : NSManagedObject

@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) NSNumber * locusId;
@property (nonatomic, retain) NSString * relationship;
@property (nonatomic, retain) NSString * sequenceId;
@property (nonatomic, retain) NSString * sequence;
@property (nonatomic, retain) NSNumber * index;
@property (nonatomic, retain) NSString * block;

@end
