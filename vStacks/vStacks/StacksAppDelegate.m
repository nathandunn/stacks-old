//
// Created by Nathan Dunn on 9/9/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
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
    NSLog(@" in this method ");
    id currentDocument = [[StacksDocumentController sharedDocumentController] currentDocument];
    NSLog(@"current doc %@",currentDocument) ;
    NSLog(@"is not nil %i",currentDocument!=nil) ;
    return [[StacksDocumentController sharedDocumentController] currentDocument]!=nil;
}

@end
