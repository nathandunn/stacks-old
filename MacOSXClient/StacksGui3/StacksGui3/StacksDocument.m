//
//  StacksDocument.m
//  StacksGui3
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StacksDocument.h"

@implementation StacksDocument

// TODO: remove these in favor of NSSet loci
@synthesize locusViews;
@synthesize loci;

- (id)init
{
    self = [super init];
    if (self) {
        // Add your subclass-specific initialization here.
    }
    return self;
}

//- (void)makeWindowControllers {
//    [super makeWindowControllers];
//}


- (NSString *)windowNibName
{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"StacksDocument";
}

- (void)windowControllerDidLoadNib:(NSWindowController *)aController
{
//    NSString *aControllerName = [anIdentifier stringByAppendingString: @"ViewController"];
//    NSString *aNibName = [anIdentifier stringByAppendingString: @"View"];
//    Class aControllerClass = NSClassFromString(aControllerName);
//    [self setCurrentController: [[aControllerClass alloc] initWithNibName: aNibName bundle: [NSBundle mainBundle]]];

    [super windowControllerDidLoadNib:aController];
    // Add any code here that needs to be executed once the windowController has loaded the document's window.


}

+ (BOOL)autosavesInPlace
{
    return YES;
}

- (id)initWithLocusView:(NSMutableDictionary*)locusViews {
    if ((self = [super init])) {
        self.locusViews = locusViews;
        self.orderedLocus = [[locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
            return [obj1 integerValue] - [obj2 integerValue];
        }];

    }
    return self ;
}

- (id)initWithLoci:(NSSet*)loci {
    if ((self = [super init])) {
        self.loci = loci ;
//        self.orderedLocus = [[locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
//            return [obj1 integerValue] - [obj2 integerValue];
//        }];

    }
    return self ;
}


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
