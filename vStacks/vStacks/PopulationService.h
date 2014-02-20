//
// Created by Nathan Dunn on 2/3/14.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>


@interface PopulationService : NSObject

+ (PopulationService*)sharedInstance;

- (NSString *)validatePopmap:(NSURL *)popmapURL;
@end
