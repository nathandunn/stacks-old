//
//  ProgressController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 9/10/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "ProgressController.h"
#import "StacksConverter.h"
#import "StacksAppDelegate.h"


@implementation ProgressController

@synthesize  stacksConverter;
@synthesize  actionMessage;
@synthesize  actionTitle;
@synthesize  cancelButton;

- (id)init {
//    NSLog(@"initializaing") ;
    self = [super init];
    if (self) {
        self = [super initWithWindowNibName:@"ProgressController"];
//        _stacksAppDelegate = [[NSApplication sharedApplication] delegate];
    }

    return self;
}



- (void)windowDidLoad
{
//    NSLog(@"trying to load window");
    [super windowDidLoad];
//    NSLog(@"window really did load");
    
    // Implement this method to handle any initialization after your window controller's window has been loaded from its nib file.
}

- (IBAction)cancelCurrentAction:(id)sender {
    NSLog(@"cancelling current action from %@",sender);
//    stacksConverter.stopProcess = true ;
    stacksConverter.stopProcess = YES ;
    actionTitle.stringValue = @"Cancelling";
    actionMessage.stringValue = @"Cancelling .. please wait.";
    [cancelButton setEnabled:NO];
//    [self close];
}


@end
