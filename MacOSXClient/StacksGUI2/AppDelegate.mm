//
//  AppDelegate.m
//  StacksGUI2
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "AppDelegate.h"
#import "MasterViewController.h"
#import "StacksLoader.h"
//#import "StacksDocument.h"
//#import "sql_utilities.h"


@implementation AppDelegate



- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
    // Insert code here to initialize your application
    // 1. Create the master view controller
    self.masterViewController = [[MasterViewController alloc] initWithWindowNibName:@"MasterViewController"];
    self.loader = [[StacksLoader alloc] init];


    NSString *examplePath = @"/tmp/stacks_tut/";
    NSMutableArray *geneDocs = [self.loader loadLoci:examplePath];
//    NSMutableArray *geneDocs;


    self.masterViewController.data = geneDocs;
    [self.masterViewController showWindow:self];

    // 2. Add the view controller to the Window's content view
//    [self.contentView addSubview:self.masterViewController.view];
//    masterViewController.view.frame = ((NSView*)self.window.contentView).bounds;
    
}


@end
