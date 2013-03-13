//
//  MasterViewController.h
//  StacksGUI2
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@class StacksDocument;
@class StacksView;
@class LocusView;

@interface MasterViewController : NSWindowController <NSSplitViewDelegate>{
    IBOutlet NSSplitView *verticalSplitView ;
    IBOutlet NSSplitView *horizontalSplitView ;
    IBOutlet NSView *dividerHandleView;

}

//@property (strong) NSMutableDictionary *stacksDocuments;
@property (atomic,strong) StacksDocument* stacksDocument;
@property (atomic,strong) LocusView* selectedLocusView;
@property(atomic, strong) StacksView *selectedStacks;

@end
