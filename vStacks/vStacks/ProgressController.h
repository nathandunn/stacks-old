//
//  ProgressController.h
//  StacksGui3
//
//  Created by Nathan Dunn on 9/10/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@class StacksAppDelegate;
@class StacksConverter;

@interface ProgressController : NSWindowController{


}

@property(strong) IBOutlet NSProgressIndicator *loadProgress;
@property(weak) IBOutlet NSButton *cancelButton;
@property(weak) IBOutlet NSTextField *actionTitle ;
@property(weak) IBOutlet NSTextField *actionMessage;
- (IBAction) cancelCurrentAction: (id) sender ;

@property(weak) StacksAppDelegate *stacksAppDelegate ;

@property(nonatomic, strong) StacksConverter *stacksConverter;
@end
