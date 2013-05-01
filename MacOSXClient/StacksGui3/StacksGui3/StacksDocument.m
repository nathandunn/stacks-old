//
//  StacksDocument.m
//  StacksGui3
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StacksDocument.h"

@implementation StacksDocument

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

@end
