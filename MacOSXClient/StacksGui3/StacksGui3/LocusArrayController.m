//
//  LocusArrayController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/17/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "LocusArrayController.h"

@implementation LocusArrayController

@synthesize minSnpValue;

- (void)awakeFromNib {
    [super awakeFromNib];
    [self setSortDescriptors:[NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"locusId" ascending:YES selector:@selector(compare:)]]];

}


- (NSArray *)arrangeObjects:(NSArray *)objects {
    return [super arrangeObjects:objects];

    if (minSnpValue == 0) {
        return [super arrangeObjects:objects];
    }
    NSMutableArray *filteredObjects = [NSMutableArray arrayWithCapacity:[objects count]];

    // these are all LocusMO objects

    NSEnumerator *objectsEnumerator = [objects objectEnumerator];
    id item;

    while (item = [objectsEnumerator nextObject]) {
//        if ([[item valueForKeyPath:@"title"] rangeOfString:searchString options:NSAnchoredSearch].location != NSNotFound) {
            [filteredObjects addObject:item];
//        }
    }
    return [super arrangeObjects:filteredObjects];
}

- (IBAction)setMinSnpValue:(id)sender {
    NSTextField *value = sender;
    minSnpValue = value.intValue;
    NSLog(@"setting value %ld", minSnpValue);
}

@end
