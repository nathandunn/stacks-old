//
//  LocusView.m
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import "LocusView.h"

@implementation LocusView
@synthesize depth;
@synthesize genotypes;


- (id)initWithId:(NSInteger)locusId {
    if ((self = [super init])) {
        self.locusId = locusId;
    }
    return self;
}

// TODO: fix parental functions
// to identify the # of parents: look at "identify_parents" in genotypes.cc .
- (NSInteger)matchingParents {
    NSInteger count = 0;
//    if (_male) {
//        ++count;
//    }
//    if (_female) {
//        ++count;
//    }
    return count;
}


@end
