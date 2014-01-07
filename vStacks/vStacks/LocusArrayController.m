//
//  LocusArrayController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/17/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "LocusArrayController.h"
#import "LocusMO.h"
#import "LocusRepository.h"

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
    minBasePairs = 0;
    maxBasePairs = 1000000 * 10000;

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

    for (id item in objects) {
        if ([item isKindOfClass:[LocusMO class]]) {
            LocusMO *locusMO = (LocusMO *) item;
            NSError *error;

            if (locusMO.snpData != nil) {
                NSDictionary *snpJson = [NSJSONSerialization JSONObjectWithData:locusMO.snpData options:kNilOptions error:&error];
                NSUInteger snpCount = snpJson.count;
                if (snpCount >= minSnpValue && snpCount <= maxSnpValue) {
                    if ([locusMO.type isEqualToString:@"Population"]) {
                        if (([locusMO.basePairs integerValue] >= minBasePairs)
                                && ([locusMO.basePairs integerValue] <= maxBasePairs)
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
    NSPopUpButton *value = sender;
    if (value.stringValue.length > 0) {
        minSnpValue = value.titleOfSelectedItem.integerValue;
        NSLog(@"setting MIN value %ld", minSnpValue);
        [self rearrangeObjects];
    }
}

- (IBAction)writeMaxSnpValue:(id)sender {
    NSPopUpButton *value = sender;
    if (value.stringValue.length > 0) {
        maxSnpValue = value.titleOfSelectedItem.integerValue;
        NSLog(@"setting MAX value %ld", maxSnpValue);
        [self rearrangeObjects];
    }
}

- (IBAction)writeLocationValue:(id)sender {
    NSPopUpButton *value = sender;
    if (value.stringValue.length > 0) {
        chromosomeLocation = value.titleOfSelectedItem;
        if ([chromosomeLocation isEqualToString:@"All Chromosomes"]) {
            chromosomeLocation = nil ;
        }
//        NSLog(@"setting Location Value %@", chromosomeLocation);
        [self rearrangeObjects];
    }
}

- (IBAction)writeMinBasePairs:(id)sender {
    NSTextField *value = sender;
    if (value.stringValue.length > 0) {
        minBasePairs = value.doubleValue * 1000000;
        NSLog(@"min baise pairs %f", minBasePairs);
        [self rearrangeObjects];
    }

}

- (IBAction)writeMaxBasePairs:(id)sender {
    NSTextField *value = sender;
    if (value.stringValue.length > 0) {
        maxBasePairs = value.doubleValue * 1000000;
        NSLog(@"max base pairs %f", maxBasePairs);
        [self rearrangeObjects];
    }

}


@end
