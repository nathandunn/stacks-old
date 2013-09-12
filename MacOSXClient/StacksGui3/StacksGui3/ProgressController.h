//
//  ProgressController.h
//  StacksGui3
//
//  Created by Nathan Dunn on 9/10/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@class StacksAppDelegate;

@interface ProgressController : NSWindowController{


}

//@property(weak) IBOutlet NSPanel *progressPanel ;
//@property(weak) IBOutlet NSProgressIndicator *loadProgress;
//@property(weak) IBOutlet NSButton *cancelButton;

@property(strong) IBOutlet NSProgressIndicator *loadProgress;
@property(weak) IBOutlet NSButton *cancelButton;
@property(weak) IBOutlet NSTextField *actionTitle ;
@property(weak) IBOutlet NSTextField *actionMessage;
//@property(weak) IBOutlet NSTe*cancelButton;
- (IBAction) cancelCurrentAction: (id) sender ;

@property(weak) StacksAppDelegate *stacksAppDelegate ;

@end
