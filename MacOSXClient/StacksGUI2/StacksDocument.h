//
//  StacksDocument.h
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>

@class LocusView;

@interface StacksDocument : NSObject


@property (strong) NSMutableDictionary *locusViews;

@property(strong) NSString *path;

@property(atomic, strong) NSMutableDictionary *populationLookup;

@property(atomic, strong) NSArray *orderedLocus;

- (id)initWithLocusView:(NSMutableDictionary*)locusViews;


- (NSMutableArray *)findPopulations;
@end
