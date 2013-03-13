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
//#import "stacks.h"


@interface MasterViewController ()

@property(weak) IBOutlet NSTableView *filesTableView;
@property(weak) IBOutlet NSTableView *genotypeTableView;
@property(weak) IBOutlet NSTableView *stacksTableView;
@property(weak) IBOutlet NSTextField *locusDetail;
@property(weak) IBOutlet NSTextField *consensusDetail;

@end

@implementation MasterViewController

@synthesize stacksDocument;

// -------------------------------------------------------------------------------
//	awakeFromNib:
// -------------------------------------------------------------------------------
- (void)awakeFromNib {
    [verticalSplitView setDelegate:self];    // we want a chance to affect the vertical split view coverage
    [_genotypeTableView setTarget:self];
    [_genotypeTableView setAllowsColumnSelection:TRUE];
    [_genotypeTableView setSelectionHighlightStyle:NSTableViewSelectionHighlightStyleNone];
    [_genotypeTableView setAction:@selector(genotypeSelected:)];
}


// -------------------------------------------------------------------------------
//	splitView:effectiveRect:effectiveRect:forDrawnRect:ofDividerAtIndex
// -------------------------------------------------------------------------------
- (NSRect)splitView:(NSSplitView *)splitView effectiveRect:(NSRect)proposedEffectiveRect forDrawnRect:(NSRect)drawnRect ofDividerAtIndex:(NSInteger)dividerIndex {
    NSRect effectiveRect = drawnRect;

    if (splitView == verticalSplitView) {
        // don't steal as much from the scroll bar as NSSplitView normally would
        effectiveRect.origin.x -= 2.0;
        effectiveRect.size.width += 6.0;

    }

    return effectiveRect;
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
            NSInteger count = [locusView genotypes];
            NSInteger rows = count / 10;
            if (count % 10 > 0) {
                rows++;
            }
            return rows;
        }
        return 0;
    }
    else
    if ([[tableView identifier] isEqualToString:@"LocusTable"]) {
        return [self.stacksDocument.locusViews count];
    }
    else
    if ([[tableView identifier] isEqualToString:@"StacksTableView"]) {
        if(self.selectedStacks==nil){
            return 0 ;
        }
        else{
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
    else
    if ([[tableView identifier] isEqualToString:@"LocusTable"]) {
        // we want data for the row . . . .
        NSArray *sortedKeys = [[self.stacksDocument.locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
            return [obj1 integerValue] - [obj2 integerValue];
        }];
        NSString *key = [sortedKeys objectAtIndexedSubscript:row];
        LocusView *locusView = [self.stacksDocument.locusViews objectForKey:key];

        // Since this is a single-column table view, this would not be necessary.
        // But it's a good practice to do it in order by remember it when a table is multicolumn.
        if ([tableColumn.identifier isEqualToString:@"IdColumn"]) {
            cellView.textField.stringValue = locusView.locusId;
        }
        else if ([tableColumn.identifier isEqualToString:@"SnpColumn"]) {
            NSMutableArray *snps = locusView.snps;
            if ([snps count] > 0) {
                cellView.textField.stringValue = [NSString stringWithFormat:@"Yes [%ldnuc]", [snps count]];
            }
            else {
                cellView.textField.stringValue = @"None";
            }
        }
        else if ([tableColumn.identifier isEqualToString:@"ParentsColumn"]) {
            cellView.textField.integerValue = [locusView matchingParents];
        }
        else if ([tableColumn.identifier isEqualToString:@"ProgenyColumn"]) {
            NSUInteger count = [[locusView progeny] count];
            cellView.textField.stringValue = [NSString stringWithFormat:@"%ld / %ld", count, count];
        }
        else if ([tableColumn.identifier isEqualToString:@"MarkerColumn"]) {
            cellView.textField.stringValue = locusView.marker;
        }
        else if ([tableColumn.identifier isEqualToString:@"RatioColumn"]) {
            cellView.textField.stringValue = @"aa: 45 (51.7%) bb:42 (48.3%)";
        }
        else if ([tableColumn.identifier isEqualToString:@"GenotypesColumn"]) {
            cellView.textField.integerValue = [locusView genotypes];
        }

        return cellView;
    }
    else
    if ([[tableView identifier] isEqualToString:@"StacksTableView"]) {
        if(self.selectedStacks!=nil){
            StacksView *stacksView = self.selectedStacks;


            if ([tableColumn.identifier isEqualToString:@"IdColumn"]) {
                if(row>2){
                    cellView.textField.integerValue = [(StackEntry *) [stacksView.stackEntries objectAtIndex:row-3] entryId];
                }
                else{
                    cellView.textField.stringValue =@"";
                }
            }
            else
            if ([tableColumn.identifier isEqualToString:@"RelationshipColumn"]) {
                switch (row){
                    case 0:
                        cellView.textField.stringValue =@"";
                        break ;
                    case 1:
                        cellView.textField.stringValue =@"consensus";
                        break ;
                    case 2:
                        cellView.textField.stringValue =@"model";
                        break ;
                    default:
                        cellView.textField.stringValue =[(StackEntry *) [stacksView.stackEntries objectAtIndex:row-3] relationship];
                }
            }
            else
            if ([tableColumn.identifier isEqualToString:@"SequenceIdColumn"]) {
                switch (row){
                    case 0:
                    case 1:
                    case 2:
                        cellView.textField.stringValue =@"";
                        break ;
                    default:
                        cellView.textField.stringValue =[(StackEntry *) [stacksView.stackEntries objectAtIndex:row-3] sequenceId];
                        cellView.textField.alignment = NSRightTextAlignment;
                }
            }
            else
            if ([tableColumn.identifier isEqualToString:@"SequenceColumn"]) {
                switch (row){
                    case 0:
                        cellView.textField.stringValue =stacksView.reference.sequence;
                        break ;
                    case 1:
                        cellView.textField.stringValue =stacksView.consensus.sequence;
                        break ;
                    case 2:
                        cellView.textField.stringValue =stacksView.model.sequence;
                        break ;
                    default:
                        NSString *sequenceString = [(StackEntry *) [stacksView.stackEntries objectAtIndex:row-3] sequence];
//                        NSMutableAttributedString *string = [self decorateSnps:sequenceString snps:stacksView.snps];
                        NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:sequenceString];

                        [string beginEditing];
                        NSNumber *snpIndex ;
                        NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                                [NSColor blueColor], NSForegroundColorAttributeName,
                                [NSColor grayColor], NSBackgroundColorAttributeName,
                                [NSFont fontWithName:@"Courier Bold" size:14.0], NSFontAttributeName,
                                nil];
                        for(snpIndex in stacksView.snps){
                            NSRange selectedRange = NSMakeRange([snpIndex intValue], 1);
                            [string setAttributes:attributes range:selectedRange];
                        }
                        [string endEditing];


                        cellView.textField.attributedStringValue = string;
                        cellView.textField.font = [NSFont fontWithName:@"Courier" size:14];
                }
            }
            else{
                NSLog(@"not sure what that column is %@",tableColumn.identifier);
                cellView.textField.stringValue = @"";
            }
        }
        return cellView;
    }
    else{
        NSLog(@"could not find table %@",[tableView identifier]);
        return cellView;
    }


}

- (NSTableCellView *)handleGenotypesTable:(NSString *)column row:(NSInteger)row cell:(NSTableCellView *)cellView {
//    NSLog(@"handling the genotypes table %@",column);
    if (self.stacksDocument != nil) {

        LocusView *locusView = self.selectedLocusView;
        NSInteger parentCount = locusView.matchingParents;
        NSInteger progenyCount = locusView.genotypes;
        NSInteger totalCount = parentCount + progenyCount;
        NSInteger totalColumnCount = 10;
        NSInteger totalRowCount = totalCount / totalColumnCount;
        NSInteger remainderColumns = totalCount % totalColumnCount;
        if (remainderColumns > 0) {
            ++totalRowCount;
        }


        // turn the row / column into an index
        NSInteger columnIndex = [[column substringFromIndex:9] integerValue];
        NSInteger progenyIndex = row * totalColumnCount + columnIndex - parentCount;

        // if a male
        if (row == 0 && [column isEqualToString:@"Genotypes1"] && [locusView hasMale]) {
            GenotypeEntry *male = locusView.male;
            cellView.textField.stringValue = [NSString stringWithFormat:@"male - %@", [male render]];
        }
        else if (row == 0 && [column isEqualToString:@"Genotypes1"] && [locusView hasFemale] && ![locusView hasMale]) {
            GenotypeEntry *female = locusView.female;
            cellView.textField.stringValue = [NSString stringWithFormat:@"female - %@", [female render]];
        }
        else if (row == 0 && [column isEqualToString:@"Genotypes2"] && [locusView hasMale] && [locusView hasFemale]) {
            GenotypeEntry *female = locusView.female;
            cellView.textField.stringValue = [NSString stringWithFormat:@"female - %@", [female render]];
        }
        else if (progenyCount > progenyIndex) {
            GenotypeEntry *genotypeEntry = [locusView.progeny objectAtIndex:progenyIndex];
            cellView.textField.stringValue = [NSString stringWithFormat:@"%ld %@", (long) genotypeEntry.entryId, [genotypeEntry render]];
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
    LocusView *locusView= [self.stacksDocument.locusViews objectForKey:key];
    return locusView;
}

- (void)handleSelectedLocus:(LocusView*) locus{

    if (locus != nil) {
        NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:locus.consensus];

        // TODOL need to store snps correctly . . . as SnpsView or as vector<SNP>:: in LocusView
//        NSMutableArray* snps = locus.snps;
//        for(int i =0 ; i < [snps count]; i++){
//           SNP *snp = [[snps objectAtIndex:i] pointerValue];
//        }

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
        [self.locusDetail setStringValue:locus.locusId];
        [self.consensusDetail setAttributedStringValue:string];

        [self.genotypeTableView reloadData];

    }
    else {
        self.stacksDocument = nil ;
        [self.locusDetail setStringValue:@"Error"];
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
//    else
//    if ([tableName isEqualToString:@"GenotypeTableView"]) {
//        NSUInteger progenyIndex = [self getSelectedGenotype];
//        self.selectedStacks = [self loadStacksForProgeny:progenyIndex];
//        [self showTags:self.selectedStacks];
//        
//        NSLog(@"genotype table selected");
//    }
//    else {
//        NSLog(@"need to handle the other case ");
//    }

}

// TODO: create the snps / stacks view
// http://genome.uoregon.edu/stacks/tag.php?db=tut_radtags&batch_id=1&sample_id=28&tag_id=277
- (void)showTagsTable:(StacksView *)view {
    [self.stacksTableView reloadData];
}

- (StacksView *)loadStacksForProgeny:(NSString*)stackKey {
    StacksLoader *loader = [[StacksLoader alloc] init];
    StacksView *stacksView = [loader loadStacksView:stackKey atPath:@"/tmp/stacks_tut"];

//    StacksView *stacksView = [[StacksView alloc] init];
    // parse the tags file based on the index
    return stacksView;
}

-(NSMutableAttributedString *)decorateSnps:(NSString *)sequenceString snps:(NSMutableArray *) snps{

}


// TODO: handle genotype selection
- (void)genotypeSelected:(id)tableView {
//    NSTableCellView *tableCellView = [tableView selectedCell];
    NSInteger rowNumber = [_genotypeTableView clickedRow];
    NSInteger columnNumber = [_genotypeTableView clickedColumn];
    if(rowNumber<0 || columnNumber<0){
        NSLog(@"invalid selection") ;
        self.selectedStacks = nil  ;
        return ;
    }
    // get the array number

    LocusView *locusView = self.selectedLocusView;
    if(rowNumber==0 && columnNumber==0 && locusView.hasMale){
        self.selectedStacks = [self loadStacksForProgeny:@"male"];
    }
    else
    if((rowNumber==0 && columnNumber==0 && !locusView.hasMale && locusView.hasFemale)
            ||
       (rowNumber==0 && columnNumber==1 && locusView.hasMale && locusView.hasFemale)
            ){
        self.selectedStacks = [self loadStacksForProgeny:@"female"];
    }
    else{
        NSUInteger totalColumnCount = 10;
        NSInteger parentCount = locusView.matchingParents;

        int index = rowNumber*totalColumnCount + columnNumber - parentCount;
        if(index+1 < [locusView genotypes]){
            GenotypeEntry *entry = (GenotypeEntry *) [locusView.progeny objectAtIndex:index+1] ;
            self.selectedStacks = [self loadStacksForProgeny:[NSString stringWithFormat:@"%ld",[entry entryId]]];
        }
        else{
            NSLog(@"invalid selection") ;
            self.selectedStacks = nil ;
        }
    }
    [self showTagsTable:self.selectedStacks];
}

@end


