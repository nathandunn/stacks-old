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

@interface MasterViewController : NSWindowController <NSSplitViewDelegate>{
    IBOutlet NSSplitView *verticalSplitView ;
    IBOutlet NSSplitView *horizontalSplitView ;
    IBOutlet NSView *dividerHandleView;

    IBOutlet NSArrayController *genotypesController ;

}

-(IBAction) genotypeSelected:(id) sender ;

@property (atomic,strong) StacksDocument* stacksDocument;
@property (atomic,strong) LocusView* selectedLocusView;
@property(atomic, strong) StacksView *selectedStacks;
@property NSInteger selectedPopulation ;



@property (strong) NSMutableArray *selectedGenotypes ;

@property(nonatomic, copy) NSString *previousStacksName;
@end
