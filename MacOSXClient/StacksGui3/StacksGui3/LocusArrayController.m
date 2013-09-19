//
//  LocusArrayController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/17/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "LocusArrayController.h"
#import "LocusMO.h"

@implementation LocusArrayController {

}


@synthesize minSnpValue;
@synthesize maxSnpValue;


- (void)awakeFromNib {
    [super awakeFromNib];
    [self setSortDescriptors:[NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"locusId" ascending:YES selector:@selector(compare:)]]];

}


- (NSArray *)arrangeObjects:(NSArray *)objects {
//    return [super arrangeObjects:objects];

    NSLog(@"calculating objects iwth minSnp %ld and maxSnp %ld for objects %ld", minSnpValue, maxSnpValue, objects.count);

    if (minSnpValue == 0 && maxSnpValue == NSIntegerMax) {
        return [super arrangeObjects:objects];
    }
    NSMutableArray *filteredObjects = [NSMutableArray arrayWithCapacity:[objects count]];

    // these are all LocusMO objects

    // TODO: objects are not all a type of locusMO

//    while (item = [objectsEnumerator nextObject]) {
    for (LocusMO *locusMO in objects) {
//        NSInteger snpCount = locusMO.snps.count;
        if (locusMO.snps != nil) {

            if (locusMO.snps.count >= minSnpValue && locusMO.snps.count <= maxSnpValue) {
                NSLog(@"snpCount %ld >= %ld && %ld <= %ld", locusMO.snps.count, minSnpValue, locusMO.snps.count, maxSnpValue);
//        if ([[item valueForKeyPath:@"title"] rangeOfString:searchString options:NSAnchoredSearch].location != NSNotFound) {
                [filteredObjects addObject:locusMO];
//        }
            }
        }
    }
    NSLog(@"filtered objects left %ld", filteredObjects.count);

    return [super arrangeObjects:filteredObjects];
}

- (IBAction)setMinSnpValue:(id)sender {
    NSTextField *value = sender;
    if (value.stringValue.length > 0) {
        minSnpValue = value.intValue;
        NSLog(@"setting MIN value %ld", minSnpValue);
        [self rearrangeObjects];
    }
}

- (IBAction)setMaxSnpValue:(id)sender {
    NSTextField *value = sender;
    if (value.stringValue.length > 0) {
        maxSnpValue = value.intValue;
        NSLog(@"setting MAX value %ld", maxSnpValue);
        [self rearrangeObjects];
    }
}

@end
