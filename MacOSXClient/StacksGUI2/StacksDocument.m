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

- (NSUInteger)childCount {
    return self.locusViews.count;
}

- (id)childAtIndex:(NSUInteger) index {
    NSArray *sortedKeys = [self.locusViews.allKeys sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
        return [obj1 integerValue] - [obj2 integerValue];
    }];
    NSString *key = [sortedKeys objectAtIndexedSubscript:index];
    return [self.locusViews objectForKey:key];
//    return nil;
}

- (BOOL)isLeaf {
    return self.childCount ==0 ;
}


@synthesize locusViews;

@end

