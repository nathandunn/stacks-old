//
// Created by Nathan Dunn on 4/1/14.
// Copyright (c) 2014 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>

@class GenericHashMO;


@interface GenericHashRepository : NSObject

+ (GenericHashRepository *)sharedInstance;

- (GenericHashMO *)store:(NSManagedObjectContext *)context key:(NSString *)key dictionary:(NSMutableDictionary *)dictionary type:(NSString *)type;

- (NSDictionary *)getDictionary:(NSManagedObjectContext *)context forKey:(NSString *)key;

- (void)detachAll;
@end