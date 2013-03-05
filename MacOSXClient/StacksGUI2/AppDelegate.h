//
//  AppDelegate.h
//  StacksGUI2
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>
//#import "MasterViewController.h"

@class MasterViewController ;
@class StacksLoader;

@interface AppDelegate : NSObject <NSApplicationDelegate>{

//    MasterViewController *masterViewController ;
}

//@property (assign) IBOutlet NSWindow *window;
@property (nonatomic,strong) IBOutlet MasterViewController *masterViewController;
@property (nonatomic,strong) StacksLoader *loader;

@end
