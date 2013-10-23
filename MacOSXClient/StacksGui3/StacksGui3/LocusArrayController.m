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
@synthesize chromosomeLocation;
@synthesize minBasePairs;
@synthesize maxBasePairs;

- (void)awakeFromNib {
    [super awakeFromNib];
    // apparently this is called before  / instead of init
    minSnpValue = 0;
    maxSnpValue = 1000;
    chromosomeLocation = nil ;
    minBasePairs = -1;
    maxBasePairs = -1;
    [self setSortDescriptors:[NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"locusId" ascending:YES selector:@selector(compare:)]]];
}

//- (NSInteger) countAll{
//    return 7 ;
//}
- (NSArray *)snpValues:(NSArray *)objects {
    NSMutableArray *filteredObjects = [NSMutableArray arrayWithCapacity:[objects count]];

    return filteredObjects;
}

- (NSArray *)allObjects:(NSArray *)objects {
    return objects;
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
    for (id item in objects) {
        if ([item isKindOfClass:[LocusMO class]]) {
            LocusMO *locusMO = (LocusMO *) item;
            if (locusMO.snps != nil) {
                if (locusMO.snps.count >= minSnpValue && locusMO.snps.count <= maxSnpValue) {
                    if ([locusMO.type isEqualToString:@"Population"]) {
                        if ((minBasePairs < 0 || [locusMO.basePairs integerValue] > minBasePairs)
                                && (maxBasePairs < 0 || [locusMO.basePairs integerValue] < maxBasePairs)
                                && (chromosomeLocation == nil  || [locusMO.chromosome isEqualToString:chromosomeLocation])
                                ) {
                            [filteredObjects addObject:locusMO];
                        }
                    }
                    else {
                        [filteredObjects addObject:locusMO];
                    }
                }
            }
        }
    }
    NSLog(@"filtered objects left %ld", filteredObjects.count);

    return [super arrangeObjects:filteredObjects];
}

- (IBAction)writeMinSnpValue:(id)sender {
    NSTextField *value = sender;
    if (value.stringValue.length > 0) {
        minSnpValue = value.intValue;
        NSLog(@"setting MIN value %ld", minSnpValue);
        [self rearrangeObjects];
    }
}

- (IBAction)writeMaxSnpValue:(id)sender {
    NSTextField *value = sender;
    if (value.stringValue.length > 0) {
        maxSnpValue = value.intValue;
        NSLog(@"setting MAX value %ld", maxSnpValue);
        [self rearrangeObjects];
    }
}

- (IBAction)writeLocationValue:(id)sender {
    NSPopUpButton *value = sender;
    if (value.stringValue.length > 0) {
        chromosomeLocation = value.titleOfSelectedItem;
        if([chromosomeLocation isEqualToString:@"All Chromosomes"]){
            chromosomeLocation= nil ;
        }
//        NSLog(@"setting Location Value %@", chromosomeLocation);
        [self rearrangeObjects];
    }
}

- (IBAction)writeMinBasePairs:(id)sender {
    NSTextField *value = sender;
    if (value.stringValue.length > 0) {
        minBasePairs = value.integerValue;
        NSLog(@"min baise pairs %ld", minBasePairs);
        [self rearrangeObjects];
    }

}

- (IBAction)writeMaxBasePairs:(id)sender {
    NSTextField *value = sender;
    if (value.stringValue.length > 0) {
        maxBasePairs = value.integerValue;
        NSLog(@"max base pairs %ld", maxBasePairs);
        [self rearrangeObjects];
    }

}


@end
