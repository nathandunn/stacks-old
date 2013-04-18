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
#import "StacksDocument.h"
#import "LocusView.h"
//#import "StacksDocument.h"
//#import "sql_utilities.h"


@implementation AppDelegate


- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
    // Insert code here to initialize your application
    // 1. Create the master view controller


    self.loader = [[StacksLoader alloc] init];

    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *examplePath = @"/tmp/stacks_tut/";
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];
    if(existsAtPath){
        [self loadApplication:examplePath];
    }
    else{
        NSLog(@"%@ does not exist.",examplePath);
    }


}

- (void)loadApplication:(NSString*)path{
    StacksDocument *stacksDocument = [self.loader loadLociAndGenotypes:path];
    self.masterViewController = [[MasterViewController alloc] initWithWindowNibName:@"MasterViewController"];
    self.masterViewController.stacksDocument = stacksDocument;
    [self.masterViewController showWindow:self];
}

-(IBAction)openDocument:(id)sender{
    int result;
    NSOpenPanel *oPanel = [NSOpenPanel openPanel];
//    NSArray *fileTypes = [NSArray arrayWithObject:@"td"];
//    [oPanel setAllowedFileTypes:fileTypes];
//    [oPanel setDirectoryURL:<#(NSURL *)url#>];
    [oPanel setAllowsMultipleSelection:NO];
    [oPanel setCanChooseDirectories:YES];
    [oPanel setCanChooseFiles:NO];

    result = [oPanel runModal];
    if (result == NSOKButton) {
        NSURL *url = [oPanel URL];
        NSString *pathToOpen = [[url path] stringByAppendingString:@"/"];
        NSLog(@"trying to poen file: %@",pathToOpen);
        [self loadApplication:pathToOpen];
    }
}


@end
