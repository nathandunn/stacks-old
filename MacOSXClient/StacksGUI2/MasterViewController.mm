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
#import "GenotypeEntry.h"
#import "StacksView.h"
#import "StacksLoader.h"
#import "StackEntry.h"
#import "LocusCell.h"
#import "stacks.h"
#import "SnpView.h"
#import "GenotypeCell.h"
//#import "stacks.h"


@interface MasterViewController ()

@property(weak) IBOutlet NSTableView *filesTableView;
@property(weak) IBOutlet NSTableView *genotypeTableView;
@property(weak) IBOutlet NSTableView *stacksTableView;
//@property(weak) IBOutlet NSTextField *locusDetail;
//@property(weak) IBOutlet NSTextField *consensusDetail;
@property(strong) IBOutlet NSWindow *mainWindow;

@end

@implementation MasterViewController

@synthesize stacksDocument;

// -------------------------------------------------------------------------------
//	awakeFromNib:
// -------------------------------------------------------------------------------
- (void)awakeFromNib {
    [verticalSplitView setDelegate:self];    // we want a chance to affect the vertical split view coverage
    [self.genotypeTableView setTarget:self];
    [self.genotypeTableView setAllowsColumnSelection:TRUE];
    [self.genotypeTableView setSelectionHighlightStyle:NSTableViewSelectionHighlightStyleNone];
    [self.genotypeTableView setAction:@selector(genotypeSelected:)];
    self.mainWindow.backgroundColor = [NSColor whiteColor];

    self.stacksLoader = [[StacksLoader alloc] init];

}


// -------------------------------------------------------------------------------
//	splitView:effectiveRect:effectiveRect:forDrawnRect:ofDividerAtIndex
// -------------------------------------------------------------------------------
- (NSRect)splitView:(NSSplitView *)splitView effectiveRect:(NSRect)proposedEffectiveRect forDrawnRect:(NSRect)drawnRect ofDividerAtIndex:(NSInteger)dividerIndex {
//    NSRect effectiveRect = drawnRect;
//
//    if (splitView == verticalSplitView) {
//        // don't steal as much from the scroll bar as NSSplitView normally would
//        effectiveRect.origin.x -= 2.0;
//        effectiveRect.size.width += 6.0;
//
//    }
//
//    return effectiveRect;
    return NSZeroRect;
}

// -------------------------------------------------------------------------------
//	splitView:additionalEffectiveRectOfDividerAtIndex:dividerIndex:
// -------------------------------------------------------------------------------
- (NSRect)splitView:(NSSplitView *)splitView additionalEffectiveRectOfDividerAtIndex:(NSInteger)dividerIndex {
    // we have a divider handle next to one of the split views in the window
    if (splitView == verticalSplitView)
        return [dividerHandleView convertRect:[dividerHandleView bounds] toView:splitView];
    else
        return NSZeroRect;
}


// -------------------------------------------------------------------------------
//	constrainMinCoordinate:proposedCoordinate:index
// -------------------------------------------------------------------------------
- (CGFloat)splitView:(NSSplitView *)splitView constrainMinCoordinate:(CGFloat)proposedCoordinate ofSubviewAt:(NSInteger)index {
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
    if ([[tableView identifier] isEqualToString:@"GenotypeTableView"]) {
        if (self.selectedLocusView != nil) {
            LocusView *locusView = self.selectedLocusView;
//            NSInteger count = [locusView genotypes];
            NSInteger count = locusView.genotypes.count;
            NSInteger rows = count / 10;
            if (count % 10 > 0) {
                rows++;
            }
            return rows;
        }
        return 0;
    }
    else if ([[tableView identifier] isEqualToString:@"LocusTable"]) {
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

- (NSView *)tableView:(NSTableView *)tableView viewForTableColumn:(NSTableColumn *)tableColumn row:(NSInteger)row {


    // Get a new ViewCell
    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];

    if ([[tableView identifier] isEqualToString:@"GenotypeTableView"]) {
//        NSLog(@"is a genotype table with column identifier %@",[tableColumn identifier]);
        return [self handleGenotypesTable:(NSString *) tableColumn.identifier row:(NSInteger) row cell:(NSTableCellView *) cellView];
    }
    else if ([[tableView identifier] isEqualToString:@"LocusTable"]) {
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
//            NSLog(@"snp index %d",snpView);
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

- (NSView *)handleStacksTable:(NSTableColumn *)tableColumn row:(NSInteger)row cell:(NSTableCellView *)cellView {
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
                    cellView.textField.stringValue = stacksView.reference.sequence;

                }
                    break;
                case 1: {
//                    cellView.textField.stringValue = stacksView.consensus.sequence;
                    NSString *consensusString = stacksView.consensus.sequence;
                    cellView.textField.attributedStringValue = [self createSnpsView:consensusString snps:stacksView.snps];
                    cellView.textField.font = [NSFont fontWithName:@"Courier" size:14];
                }
                    break;
                case 2: {
                    cellView.textField.stringValue = stacksView.model.sequence;
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

- (NSAttributedString *)createSnpsView:(NSString *)sequenceString snps:(NSMutableArray *)snps {
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

- (NSTableCellView *)handleGenotypesTable:(NSString *)column row:(NSInteger)row cell:(GenotypeCell *)cellView {
//    NSLog(@"handling the genotypes table %@",column);
    if (self.stacksDocument != nil) {

        LocusView *locusView = self.selectedLocusView;
        NSInteger progenyCount = locusView.genotypes.count;

        // TODO: this should come from the file system, nowhere else
// to identify the # of parents: look at "identify_parents" in genotypes.cc .
        NSInteger totalColumnCount = 10;
//        NSInteger remainderColumns = totalCount % totalColumnCount;


        // turn the row / column into an index
        NSInteger columnIndex = [[column substringFromIndex:9] integerValue];
        NSInteger progenyIndex = row * totalColumnCount + columnIndex;

        if (progenyCount > progenyIndex) {
            NSArray *sortedKeys = [[locusView.genotypes allKeys] sortedArrayUsingComparator:^(NSString *obj1, NSString *obj2) {
                return [obj1 compare:obj2];
            }];
            NSString *key = [sortedKeys objectAtIndex:progenyIndex - 1];
            GenotypeEntry *genotypeEntry = [locusView.genotypes valueForKey:key];
//            cellView.textField.stringValue = [NSString stringWithFormat:@"%@  %@", genotypeEntry.name, [genotypeEntry render]];

            NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:genotypeEntry.name];


            NSMutableAttributedString *haplotypes = [[NSMutableAttributedString alloc] initWithString:@""];

            [haplotypes beginEditing];

            for (NSString *haplotype in genotypeEntry.haplotypes) {
                NSMutableAttributedString *hapString = [[NSMutableAttributedString alloc] initWithString:haplotype];
                [haplotypes appendAttributedString:hapString];
            }

            [haplotypes endEditing];

            cellView.textField.attributedStringValue = string;
//            cellView.haplotypes.attributedStringValue = haplotypes;
        }
        else {
            cellView.textField.stringValue = @"";
        }
    }
    else {
        cellView.textField.stringValue = @"";
    }

    return cellView;
}

- (LocusView *)findSelectedLocus {
    NSInteger selectedRow = [self.filesTableView selectedRow];
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
//        [self.locusDetail setStringValue:locus.locusId];
//        [self.consensusDetail setAttributedStringValue:string];

        [self.genotypeTableView reloadData];

    }
    else {
        self.stacksDocument = nil ;
//        [self.locusDetail setStringValue:@"Error"];
    }


    // update the genotypes table
    [self.genotypeTableView reloadData];


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


// TODO: handle genotype selection
- (void)genotypeSelected:(id)tableView {
    NSInteger rowNumber = [_genotypeTableView clickedRow];
    NSInteger columnNumber = [_genotypeTableView clickedColumn];
    if (rowNumber < 0 || columnNumber < 0) {
        NSLog(@"invalid selection");
        self.selectedStacks = nil ;
        return;
    }
    // get the array number

    LocusView *locusView = self.selectedLocusView;
    NSUInteger totalColumnCount = 10;

    int index = rowNumber * totalColumnCount + columnNumber;
    NSLog(@"loading genotypes %d", locusView.genotypes.count);

    if (index + 1 < locusView.genotypes.count) {
        NSString *key = [[locusView.genotypes allKeys] objectAtIndex:index + 1];
        GenotypeEntry *genotypeEntry = [locusView.genotypes valueForKey:key];
        NSLog(@"entry ID: %d", genotypeEntry.sampleId);

        NSLog(@"loading %@ tag - %d", genotypeEntry.name, genotypeEntry.tagId);
        StacksView *stacksView = [_stacksLoader loadStacksView:genotypeEntry.name atPath:@"/tmp/stacks_tut/" forTag:genotypeEntry.tagId locus:locusView];
        self.selectedStacks = stacksView;
    }
    else {
        NSLog(@"invalid selection");
        self.selectedStacks = nil ;
    }
    [self showTagsTable:self.selectedStacks];
}

@end


