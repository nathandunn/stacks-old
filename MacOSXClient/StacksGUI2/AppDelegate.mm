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


- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
    // Insert code here to initialize your application
    // 1. Create the master view controller
    self.masterViewController = [[MasterViewController alloc] initWithWindowNibName:@"MasterViewController"];
    self.loader = [[StacksLoader alloc] init];


    NSString *examplePath = @"/tmp/stacks_tut/";
    NSMutableDictionary *geneDocs = [self.loader loadLoci:examplePath];
//    NSMutableDictionary *genotypes = [self.loader loadGenotypes:examplePath withLoci:geneDocs];

    // on each call of the genotypes . . . the similarly linked genotype will also be linked
//    NSLog(@"number of genotypes %d",[genotypes count]);

//        load
//    NSMutableArray *geneDocs;


    self.masterViewController.stacksDocuments = geneDocs;
    [self.masterViewController showWindow:self];

    // 2. Add the view controller to the Window's content view
//    [self.contentView addSubview:self.masterViewController.view];
//    masterViewController.view.frame = ((NSView*)self.window.contentView).bounds;

}


@end
