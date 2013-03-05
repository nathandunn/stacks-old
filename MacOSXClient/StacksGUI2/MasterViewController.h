//
//  MasterViewController.h
//  StacksGUI2
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@class StacksDocument;

@interface MasterViewController : NSWindowController <NSSplitViewDelegate>{
    IBOutlet NSSplitView *verticalSplitView ;
    IBOutlet NSSplitView *horizontalSplitView ;
    IBOutlet NSView *dividerHandleView;
}

@property (strong) NSMutableDictionary *stacksDocuments;
@property (weak) StacksDocument* selectedStacksDocument;


@end
