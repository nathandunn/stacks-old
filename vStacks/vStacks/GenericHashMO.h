//
//  GenericHash.h
//  vStacks
//
//  Created by Nathan Dunn on 3/26/14.
//  Copyright (c) 2014 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>


@interface GenericHashMO : NSManagedObject

@property(nonatomic, retain) NSData *dataValue;
@property(nonatomic, retain) NSString *key;
@property(nonatomic, retain) NSString *stringValue;
@property(nonatomic, retain) NSString *type;

@end
