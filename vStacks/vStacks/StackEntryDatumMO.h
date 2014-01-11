//
//  StackEntryDatumMO.h
//  vStacks
//
//  Created by Nathan Dunn on 1/11/14.
//  Copyright (c) 2014 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>


@interface StackEntryDatumMO : NSManagedObject

@property (nonatomic, copy) NSString * name;
@property (nonatomic, copy) NSNumber * sampleId;
@property (nonatomic, copy) NSData * stackData;
@property (nonatomic, copy) NSNumber * tagId;

@end
