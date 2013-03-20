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
@class StacksLoader;
@class GenotypeEntry;

@interface MasterViewController : NSWindowController <NSSplitViewDelegate>{
}

//@property (strong) NSMutableDictionary *stacksDocuments;
@property (atomic,strong) StacksDocument* stacksDocument;
@property (atomic,strong) LocusView* selectedLocusView;
@property (atomic,strong) GenotypeEntry* selectedGenotype;
@property(atomic, strong) StacksView *selectedStacks;
@property(atomic, strong) StacksLoader *stacksLoader;

@end
