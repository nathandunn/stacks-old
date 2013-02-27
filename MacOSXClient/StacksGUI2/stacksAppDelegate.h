//
//  stacksAppDelegate.h
//  StacksGUI2
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "MasterViewController.h"

@interface stacksAppDelegate : NSObject <NSApplicationDelegate>

@property (assign) IBOutlet NSWindow *window;
@property (nonatomic,strong) IBOutlet MasterViewController *masterViewController;

@end
