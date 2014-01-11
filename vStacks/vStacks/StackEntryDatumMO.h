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

@property (nonatomic, retain) NSString * name;
@property (nonatomic, retain) NSNumber * sampleId;
@property (nonatomic, retain) id stackData;
@property (nonatomic, retain) NSNumber * tagId;

@end
