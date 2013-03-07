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
//#import "stacks.h"


@interface MasterViewController ()

@property(weak) IBOutlet NSTableView *filesTableView;
@property(weak) IBOutlet NSTableView *genotypeTableView;
@property(weak) IBOutlet NSTextField *locusDetail;
@property(weak) IBOutlet NSTextField *consensusDetail;

- (IBAction)selectGenotype:(NSTextFieldCell *)sender;

@end

@implementation MasterViewController
@synthesize selectedGenotype = _selectedGenotype;


// -------------------------------------------------------------------------------
//	awakeFromNib:
// -------------------------------------------------------------------------------
- (void)awakeFromNib {
    [verticalSplitView setDelegate:self];    // we want a chance to affect the vertical split view coverage
    [_genotypeTableView setTarget:self];
    [_genotypeTableView setAllowsColumnSelection:TRUE];
    [_genotypeTableView setDoubleAction:@selector(doubleClick:)];
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
        if (self.selectedStacksDocument != nil) {
            LocusView *locusView = self.selectedStacksDocument.locusData;
            NSInteger count = [locusView genotypes];
            NSInteger rows = count / 10;
            if (count % 10 > 0) {
                rows++;
            }
            return rows;

        }
        return 0;
    }
    else {
        return [self.stacksDocuments count];
    }
}

- (NSView *)tableView:(NSTableView *)tableView viewForTableColumn:(NSTableColumn *)tableColumn row:(NSInteger)row {


    // Get a new ViewCell
    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];

    if ([[tableView identifier] isEqualToString:@"GenotypeTableView"]) {
//        NSLog(@"is a genotype table with column identifier %@",[tableColumn identifier]);
        return [self handleGenotypesTable:(NSString *) tableColumn.identifier row:(NSInteger) row cell:(NSTableCellView *) cellView];
    }
    else {
        // we want data for the row . . . .
        NSArray *sortedKeys = [[self.stacksDocuments allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
            return [obj1 integerValue] - [obj2 integerValue];
        }];
        NSString *key = [sortedKeys objectAtIndexedSubscript:row];
        StacksDocument *stacksDoc = [self.stacksDocuments objectForKey:key];

        // Since this is a single-column table view, this would not be necessary.
        // But it's a good practice to do it in order by remember it when a table is multicolumn.
        if ([tableColumn.identifier isEqualToString:@"IdColumn"]) {
            cellView.textField.stringValue = stacksDoc.locusData.locusId;
        }
        else if ([tableColumn.identifier isEqualToString:@"SnpColumn"]) {
            NSMutableArray *snps = stacksDoc.locusData.snps;
            if ([snps count] > 0) {
                cellView.textField.stringValue = [NSString stringWithFormat:@"Yes [%ldnuc]", [snps count]];
            }
            else {
                cellView.textField.stringValue = @"None";
            }
        }
        else if ([tableColumn.identifier isEqualToString:@"ParentsColumn"]) {
            cellView.textField.integerValue = [stacksDoc.locusData matchingParents];
        }
        else if ([tableColumn.identifier isEqualToString:@"ProgenyColumn"]) {
            NSUInteger count = [[stacksDoc.locusData progeny] count];
            cellView.textField.stringValue = [NSString stringWithFormat:@"%ld / %ld", count, count];
        }
        else if ([tableColumn.identifier isEqualToString:@"MarkerColumn"]) {
            cellView.textField.stringValue = stacksDoc.locusData.marker;
        }
        else if ([tableColumn.identifier isEqualToString:@"RatioColumn"]) {
            cellView.textField.stringValue = @"aa: 45 (51.7%) bb:42 (48.3%)";
        }
        else if ([tableColumn.identifier isEqualToString:@"GenotypesColumn"]) {
            cellView.textField.integerValue = [stacksDoc.locusData genotypes];
        }

        return cellView;
    }


}

- (NSTableCellView *)handleGenotypesTable:(NSString *)column row:(NSInteger)row cell:(NSTableCellView *)cellView {
//    NSLog(@"handling the genotypes table %@",column);
    if (self.selectedStacksDocument != nil) {

        LocusView *locusView = self.selectedStacksDocument.locusData;
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

- (StacksDocument *)selectedDoc {

    NSInteger selectedRow = [self.filesTableView selectedRow];
    NSArray *sortedKeys = [[self.stacksDocuments allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
        return [obj1 integerValue] - [obj2 integerValue];
    }];
    NSString *key = [sortedKeys objectAtIndexedSubscript:selectedRow];
    StacksDocument *stacksDocument = [self.stacksDocuments objectForKey:key];
    return stacksDocument;
}

- (void)setDetailInfo:(StacksDocument *)doc {

    if (doc != nil) {
        LocusView *locus = doc.locusData;
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
        self.selectedStacksDocument = nil ;
        [self.locusDetail setStringValue:@"Error"];
    }


    // update the genotypes table
    [self.genotypeTableView reloadData];


}

- (void)tableViewSelectionDidChange:(NSNotification *)aNotification {
    NSString *tableName = [[aNotification object] identifier];
    if ([tableName isEqualToString:@"LocusTable"]) {
        self.selectedStacksDocument = [self selectedDoc];
        // Update info
        [self setDetailInfo:self.selectedStacksDocument];
    }
//    else
//    if ([tableName isEqualToString:@"GenotypeTableView"]) {
//        NSUInteger progenyIndex = [self getSelectedGenotype];
//        self.selectedGenotype = [self loadStacksForProgeny:progenyIndex];
//        [self showTags:self.selectedGenotype];
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

}
// TODO: somehow get the selected path
- (void)tableView: (NSTableView *)tableView didSelectRowAtIndexPath: (NSIndexPath *)indexPath {
    NSLog(@"getting the path . . . %@",indexPath);
//    [tableView cel]
//    NSTableCellView *cell = [tableView cellForRowAtIndexPath:indexPath];
//    NSTableCellView *cell = [tableView ];

//    NSString *newText = [array objectAtIndex:row];
//    textbox.text = newtext;
}

- (StacksView *)loadStacksForProgeny:(NSUInteger)i {
    StacksView *stacksView = [[StacksView alloc] init];
    // parse the tags file based on the index
    return stacksView;
}


// TODO: handle genotype selection
- (void)doubleClick:(id)doubleClick {
    NSInteger rowNumber = [_genotypeTableView clickedRow];
    NSInteger columnNumber = [_genotypeTableView clickedColumn];
    NSLog(@"clicked row %ld, column %ld", rowNumber, columnNumber);
    NSTableCellView *tableCellView = [_genotypeTableView selectedCell];
//    NSCell *tableCellView = [_genotypeTableView preparedCellAtColumn:columnNumber row:rowNumber];
    NSLog(@"value of cell %@", [[tableCellView textField] stringValue]);


    NSUInteger progenyIndex = 3 ;
    self.selectedGenotype = [self loadStacksForProgeny:progenyIndex];
    [self showTagsTable:self.selectedGenotype];
}

@end


