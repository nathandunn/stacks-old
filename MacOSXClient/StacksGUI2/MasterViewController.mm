//
//  MasterViewController.m
//  StacksGUI2
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "MasterViewController.h"
#import "StacksDocument.h"
#import "LocusView.h"
#import "StacksView.h"
#import "StacksLoader.h"
#import "StackEntry.h"
#import "LocusCell.h"
#import "SnpView.h"
#import "GenotypeEntry.h"


@interface MasterViewController ()

@property(weak) IBOutlet NSTableView *stacksTableView;
@property(weak) IBOutlet NSTableView *locusTableView;
@property(strong) IBOutlet NSWindow *mainWindow;
@property(strong) StacksLoader *stacksLoader;

@end

@implementation MasterViewController

@synthesize stacksDocument;
@synthesize selectedGenotypes;

// -------------------------------------------------------------------------------
//	awakeFromNib:
// -------------------------------------------------------------------------------
- (void)awakeFromNib {
    [verticalSplitView setDelegate:self];    // we want a chance to affect the vertical split view coverage
    self.mainWindow.backgroundColor = [NSColor whiteColor];
    self.stacksLoader = [[StacksLoader alloc] init];

    [genotypesController addObserver:self forKeyPath:@"selectionIndexes" options:NSKeyValueObservingOptionNew context:nil];
}


// -------------------------------------------------------------------------------
//	splitView:effectiveRect:effectiveRect:forDrawnRect:ofDividerAtIndex
// -------------------------------------------------------------------------------
- (NSRect)splitView:(NSSplitView *)splitView
      effectiveRect:
              (NSRect)proposedEffectiveRect
       forDrawnRect:
               (NSRect)drawnRect
   ofDividerAtIndex:
           (NSInteger)dividerIndex {
    return NSZeroRect;
}

// -------------------------------------------------------------------------------
//	splitView:additionalEffectiveRectOfDividerAtIndex:dividerIndex:
// -------------------------------------------------------------------------------
- (NSRect)                    splitView:(NSSplitView *)splitView
additionalEffectiveRectOfDividerAtIndex:
        (NSInteger)dividerIndex {
    // we have a divider handle next to one of the split views in the window
    if (splitView == verticalSplitView)
        return [dividerHandleView convertRect:[dividerHandleView bounds] toView:splitView];
    else
        return NSZeroRect;
}


// -------------------------------------------------------------------------------
//	constrainMinCoordinate:proposedCoordinate:index
// -------------------------------------------------------------------------------
- (CGFloat)  splitView:(NSSplitView *)splitView
constrainMinCoordinate:
        (CGFloat)proposedCoordinate
           ofSubviewAt:
                   (NSInteger)index {
    CGFloat constrainedCoordinate = proposedCoordinate;
    if (splitView == verticalSplitView) {
        // the primary vertical split view is asking for a constrained size
        constrainedCoordinate = proposedCoordinate + 120.0;
    }
    else if (splitView == horizontalSplitView) {
        // the horizontal split view between mailboxes and activity view
        constrainedCoordinate = proposedCoordinate + 200.0;
    }

    return constrainedCoordinate;
}

- (NSInteger)numberOfRowsInTableView:(NSTableView *)tableView {
    if ([[tableView identifier] isEqualToString:@"LocusTable"]) {
        return [self.stacksDocument.locusViews count];
    }
    else if ([[tableView identifier] isEqualToString:@"StacksTableView"]) {
        if (self.selectedStacks == nil) {
            return 0;
        }
        else {
            return [self.selectedStacks rowsNeeded];
        }
    }

}

- (NSView *)tableView:(NSTableView *)tableView
   viewForTableColumn:
           (NSTableColumn *)tableColumn
                  row:
                          (NSInteger)row {

    // Get a new ViewCell
    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];

    if ([[tableView identifier] isEqualToString:@"LocusTable"]) {
        // we want data for the row . . . .
        NSArray *sortedKeys = [[self.stacksDocument.locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
            return [obj1 integerValue] - [obj2 integerValue];
        }];
        NSString *key = [sortedKeys objectAtIndexedSubscript:row];
        LocusView *locusView = [self.stacksDocument.locusViews objectForKey:key];

        LocusCell *locusCell = (LocusCell *) cellView;
        locusCell.locusId.stringValue = locusView.locusId;
        locusCell.propertyField.stringValue = [NSString stringWithFormat:@"Parents %d Progeny %d \nSNPS %d"
                , 0, locusView.genotypes.count, locusView.snps.count];


        // START FANCY
        NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:locusView.consensus];

        [string beginEditing];
        NSNumber *snpIndex;
        NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blueColor], NSForegroundColorAttributeName,
                [NSColor grayColor], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier Bold" size:14.0], NSFontAttributeName,
                nil];
        for (SnpView *snpView in locusView.snps) {
            NSRange selectedRange = NSMakeRange(snpView.column, 1);
            [string setAttributes:attributes range:selectedRange];
        }
        [string endEditing];


        locusCell.consensusField.attributedStringValue = string;

        return cellView;
    }
    else if ([[tableView identifier] isEqualToString:@"StacksTableView"]) {
        return [self handleStacksTable:(NSTableColumn *) tableColumn row:(NSInteger) row cell:(NSTableCellView *) cellView];
    }
    else {
        NSLog(@"could not find table %@", [tableView identifier]);
        return cellView;
    }


}

- (NSView *)handleStacksTable:(NSTableColumn *)tableColumn
                          row:
                                  (NSInteger)row
                         cell:
                                 (NSTableCellView *)cellView {
    if (self.selectedStacks != nil) {
        StacksView *stacksView = self.selectedStacks;


        if ([tableColumn.identifier isEqualToString:@"IdColumn"]) {
            if (row > 2) {
                cellView.textField.integerValue = [(StackEntry *) [stacksView.stackEntries objectAtIndex:row - 3] entryId];
            }
            else {
                cellView.textField.stringValue = @"";
            }
        }
        else if ([tableColumn.identifier isEqualToString:@"RelationshipColumn"]) {
            switch (row) {
                case 0:
                    cellView.textField.stringValue = @"";
                    break;
                case 1:
                    cellView.textField.stringValue = @"consensus";
                    break;
                case 2:
                    cellView.textField.stringValue = @"model";
                    break;
                default:
                    cellView.textField.stringValue = [(StackEntry *) [stacksView.stackEntries objectAtIndex:row - 3] relationship];
            }
        }
        else if ([tableColumn.identifier isEqualToString:@"SequenceIdColumn"]) {
            switch (row) {
                case 0:
                case 1:
                case 2:
                    cellView.textField.stringValue = @"";
                    break;
                default:
                    cellView.textField.stringValue = [(StackEntry *) [stacksView.stackEntries objectAtIndex:row - 3] sequenceId];
                    cellView.textField.alignment = NSRightTextAlignment;
            }
        }
        else if ([tableColumn.identifier isEqualToString:@"SequenceColumn"]) {
            switch (row) {
                case 0: {
//                    cellView.textField.stringValue = stacksView.reference.sequence;
                    cellView.textField.attributedStringValue = [self createReferenceView:stacksView.consensus.sequence.length];
                    cellView.textField.font = [NSFont fontWithName:@"Courier" size:14];
                }
                    break;
                case 1: {
                    NSString *consensusString = stacksView.consensus.sequence;
                    cellView.textField.attributedStringValue = [self createSnpsView:consensusString snps:stacksView.snps];
                    cellView.textField.font = [NSFont fontWithName:@"Courier" size:14];
                }
                    break;
                case 2: {
                    cellView.textField.stringValue = stacksView.model.sequence;
                    cellView.textField.font = [NSFont fontWithName:@"Courier" size:14];
                }
                    break;
                default: {
                    NSString *sequenceString = [(StackEntry *) [stacksView.stackEntries objectAtIndex:row - 3] sequence];
                    cellView.textField.attributedStringValue = [self createSnpsView:sequenceString snps:stacksView.snps];
                    cellView.textField.font = [NSFont fontWithName:@"Courier" size:14];
                }

            }
        }
        else {
            NSLog(@"not sure what that column is %@", tableColumn.identifier);
            cellView.textField.stringValue = @"";
        }
    }
    return cellView;
}

- (NSAttributedString *)createReferenceView:(NSUInteger)sequenceSize {
    // create a string from 0-9 for sequenceSize
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:@"123123123"];

    return string;
}

- (NSAttributedString *)createSnpsView:(NSString *)sequenceString
                                  snps:
                                          (NSMutableArray *)snps {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:sequenceString];

    [string beginEditing];
    NSNumber *snpIndex;
    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blueColor], NSForegroundColorAttributeName,
            [NSColor grayColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier Bold" size:14.0], NSFontAttributeName,
            nil];
    for (snpIndex in snps) {
        NSRange selectedRange = NSMakeRange([snpIndex intValue], 1);
        [string setAttributes:attributes range:selectedRange];
    }
    [string endEditing];
    return string;
}


- (LocusView *)findSelectedLocus {
    NSInteger selectedRow = [self.locusTableView selectedRow];
    NSArray *sortedKeys = [[self.stacksDocument.locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
        return [obj1 integerValue] - [obj2 integerValue];
    }];
    NSString *key = [sortedKeys objectAtIndexedSubscript:selectedRow];
    LocusView *locusView = [self.stacksDocument.locusViews objectForKey:key];
    return locusView;
}

- (void)handleSelectedLocus:(LocusView *)locus {

    if (locus != nil) {
        NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:locus.consensus];

        NSRange selectedRange = NSMakeRange(12, 1);

        [string beginEditing];

        [NSFont boldSystemFontOfSize:10.0];
        NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blueColor], NSForegroundColorAttributeName,
                [NSColor grayColor], NSBackgroundColorAttributeName,
                [NSFont boldSystemFontOfSize:14.0], NSFontAttributeName,
                nil];
        [string setAttributes:attributes range:selectedRange];

        [string endEditing];


        self.selectedGenotypes = [locus.genotypes allValues];
    }
    else {
        self.stacksDocument = nil ;
    }


}

- (void)tableViewSelectionDidChange:(NSNotification *)aNotification {
    NSString *tableName = [[aNotification object] identifier];
    if ([tableName isEqualToString:@"LocusTable"]) {
        self.selectedLocusView = [self findSelectedLocus];
        // Update info
        [self handleSelectedLocus:self.selectedLocusView];
    }

}

// TODO: create the snps / stacks view
// http://genome.uoregon.edu/stacks/tag.php?db=tut_radtags&batch_id=1&sample_id=28&tag_id=277
- (void)showTagsTable:(StacksView *)view {
    [self.stacksTableView reloadData];
}

- (void)observeValueForKeyPath:(NSString *)keyPath
                      ofObject:(id)object
                        change:(NSDictionary *)change
                       context:(void *)context {
    if ([keyPath isEqualTo:@"selectionIndexes"]) {
        if ([[genotypesController selectedObjects] count] > 0) {
            if ([[genotypesController selectedObjects] count] == 1) {
                GenotypeEntry *genotypeEntry = (GenotypeEntry *) [[genotypesController selectedObjects] objectAtIndex:0];
                
                if([genotypeEntry.name isEqualToString:self.previousStacksName]) {
                    return ; 
                }

                self.previousStacksName = genotypeEntry.name ;
                
//                NSLog(@"selected genotype %@ and tagID %ld", genotypeEntry.name,genotypeEntry.tagId);
                LocusView *locusView = self.selectedLocusView;

//                NSLog(@"locusView: %@",locusView.locusId);
                // TODO: should use the current path of the
                StacksView *stacksView = [self.stacksLoader loadStacksView:genotypeEntry.name atPath:@"/tmp/stacks_tut/" forTag:genotypeEntry.tagId locus:locusView];
                self.selectedStacks = stacksView;
//                NSLog(@"stacks view %@",stacksView);

//                [self showTagsTable:self.selectedStacks];
            }
        }
        else {
            self.selectedStacks = nil ;
            NSLog(@"none selected");
        }
        [self showTagsTable:self.selectedStacks];
        [self.stacksTableView reloadData];
    }
}



@end


