//
//  AppDelegate.m
//  StacksGUI2
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "AppDelegate.h"
#import "MasterViewController.h"
#import "BrowserViewController.h"
#import "StacksLoader.h"
#import "StacksDocument.h"
#import "LocusView.h"
//#import "StacksDocument.h"
//#import "sql_utilities.h"


@implementation AppDelegate


- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
    // Insert code here to initialize your application
    // 1. Create the master view controller
    self.masterViewController = [[MasterViewController alloc] initWithWindowNibName:@"MasterViewController"];
    self.browserViewController = [[BrowserViewController alloc] init];
    self.loader = [[StacksLoader alloc] init];


    NSString *examplePath = @"/tmp/stacks_tut/";
//    StacksDocument *stacksDocument = [self.loader loadLoci:examplePath];
    StacksDocument *stacksDocument = [self.loader loadLociAndGenotypes:examplePath];

    self.masterViewController.stacksDocument = stacksDocument;

    [self.masterViewController showWindow:self];
//    [self.browserViewController ];

}


@end
