//
//  AlleleMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>


@interface AlleleMO : NSManagedObject

@property (nonatomic, retain) NSNumber * allele;
@property (nonatomic, retain) NSNumber * depth;
@property (nonatomic, retain) NSNumber * ratio;

@end
