//
//  SnpMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/8/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>


@interface SnpMO : NSManagedObject

@property (nonatomic, retain) NSNumber * column;
@property (nonatomic, retain) NSNumber * lratio;
@property (nonatomic, retain) NSNumber * rank1;
@property (nonatomic, retain) NSNumber * rank2;
@property (nonatomic, retain) NSNumber * rank3;
@property (nonatomic, retain) NSNumber * rank4;
@property (nonatomic, retain) NSNumber * type;

@end
