//
// Created by Nathan Dunn on 4/15/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class LocusView;


@interface PopulationView : NSObject


@property (strong) NSString *name;
@property (strong) LocusView *parentLocus;
@property (strong) NSMutableDictionary *genotypes;

@end