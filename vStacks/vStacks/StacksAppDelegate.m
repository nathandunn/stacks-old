//
// Created by Nathan Dunn on 9/9/13.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//
//
//


#import "StacksAppDelegate.h"
#import "StacksApplicationController.h"
#import "StacksDocumentController.h"


@implementation StacksAppDelegate {


}
- (id)init {
    self = [super init];
    if (self) {
        // registers controller as the shared internal controller internally 
        StacksApplicationController *stacksApplicationController = [[StacksApplicationController alloc] init];
    }

    return self;
}


- (void)applicationWillFinishLaunching:(NSNotification *)notification {
    NSLog(@"app did finish launcing");

}

//- (BOOL)validateUserInterfaceItem:(id < NSValidatedUserInterfaceItem >)anItem{
//    NSLog(@"should be avlidating something");
//    return YES;
//}

- (BOOL) hasOpenDocument{
    id currentDocument = [[StacksDocumentController sharedDocumentController] currentDocument];
    NSLog(@"Valid Document %i",currentDocument!=nil) ;
    return [[StacksDocumentController sharedDocumentController] currentDocument]!=nil;
}

@end

