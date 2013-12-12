//
//  MatchMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 12/12/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class SampleMO;

@interface MatchMO : NSManagedObject

@property (nonatomic, retain) NSNumber * internalId;
@property (nonatomic, retain) NSNumber * externalId;
@property (nonatomic, retain) SampleMO *sample;

@end
