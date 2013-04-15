//
//  StacksDocument.m
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import "StacksDocument.h"
#import "LocusView.h"

@implementation StacksDocument

- (id)initWithLocusView:(NSMutableDictionary*)locusViews {
    if ((self = [super init])) {
        self.locusViews = locusViews;
    }
    return self ;
}

@synthesize locusViews;

- (NSMutableArray *)findPopulations {

    NSMutableArray *populations = [[NSMutableArray alloc] init];

    NSString *population ;
    for(population in self.populationLookup.allValues){
        if(![populations containsObject:population]){
            [populations addObject:population];
        }
    }

    return populations;
}

@end

