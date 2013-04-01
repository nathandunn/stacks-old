//
//  BrowserViewController.h
//  StacksGUI2
//
//  Created by Nathan Dunn on 4/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@class StacksDocument;
@class StacksLoader;

@interface BrowserViewController : NSWindowController{


}

@property (weak) NSBrowser* browser;
@property (atomic,strong) StacksDocument* stacksDocument;
@property (atomic,strong) StacksLoader *stacksLoader ;

@end
