//
//  ProgressController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 9/10/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "ProgressController.h"

@interface ProgressController ()

@end

//@synthesize  loadProgress;
//@synthesize  progressPanel;


@implementation ProgressController

- (id)init {
    NSLog(@"initializaing") ;
    self = [super init];
    if (self) {
        self = [super initWithWindowNibName:@"ProgressController"];
    }

    return self;
}


//- (id)initWithWindow:(NSWindow *)window
//{
//    self = [super initWithWindow:window];
//    if (self) {
//        // Initialization code here.
//    }
//
//    return self;
//}

- (void)windowDidLoad
{
    NSLog(@"trying to load window");
    [super windowDidLoad];
    NSLog(@"window really did load");
    
    // Implement this method to handle any initialization after your window controller's window has been loaded from its nib file.
}

- (IBAction)cancelCurrentAction:(id)sender {
    NSLog(@"cancelling current action from %@",sender);
}


@end
