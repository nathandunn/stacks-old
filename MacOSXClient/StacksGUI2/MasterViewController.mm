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
//#import "stacks.h"


@interface MasterViewController ()

@property(weak) IBOutlet NSTableView *filesTableView;
@property(weak) IBOutlet NSTableView *genotypeTableView;
@property(weak) IBOutlet NSTableView *stacksTableView;
@property(weak) IBOutlet NSTextField *locusDetail;
@property(weak) IBOutlet NSTextField *consensusDetail;

@property(weak) IBOutlet NSScrollView *snpsScrollView ;
@property(weak) IBOutlet NSScrollView *genotypesScrollView ;

@end

@implementation MasterViewController

@synthesize stacksDocument;

// -------------------------------------------------------------------------------
//	awakeFromNib:
// -------------------------------------------------------------------------------
- (void)awakeFromNib {
//    [self.snpsScrollView setHidden:TRUE];
    [self.genotypeTableView setTarget:self];
    [self.genotypeTableView setAllowsColumnSelection:TRUE];
    [self.genotypeTableView setSelectionHighlightStyle:NSTableViewSelectionHighlightStyleNone];
    [self.genotypeTableView setAction:@selector(genotypeSelected:)];

    _stacksLoader = [[StacksLoader alloc] init];
}


- (NSInteger)numberOfRowsInTableView:(NSTableView *)tableView {
    if ([[tableView identifier] isEqualToString:@"GenotypeTableView"]
            ||
            [[tableView identifier] isEqualToString:@"ProgenyLightTable"]
            ) {
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
    else if ([[tableView identifier] isEqualToString:@"PopulationTable"]) {
        NSLog(@"setting num rows for pop table");
        return 1;
    }
    else{
        NSLog(@"table not found %@",tableView.identifier);
    }

}

- (NSView *)tableView:(NSTableView *)tableView viewForTableColumn:(NSTableColumn *)tableColumn row:(NSInteger)row {


    // Get a new ViewCell
    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];

    if ([[tableView identifier] isEqualToString:@"GenotypeTableView"]) {
//        NSLog(@"is a genotype table with column identifier %@",[tableColumn identifier]);
        return [self handleGenotypesTable:(NSString *) tableColumn.identifier row:(NSInteger) row cell:(NSTableCellView *) cellView];
    }
    if ([[tableView identifier] isEqualToString:@"PopulationTable"]) {
        NSLog(@"populating the population table  %@",tableColumn.identifier);
        cellView.textField.stringValue = @"Rabbit Slough";
        return cellView;
//        return [self handleGenotypesTable:(NSString *) tableColumn.identifier row:(NSInteger) row cell:(NSTableCellView *) cellView];
    }
    else if ([[tableView identifier] isEqualToString:@"LocusTable"]) {
        // we want data for the row . . . .
        NSArray *sortedKeys = [[self.stacksDocument.locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
            return [obj1 integerValue] - [obj2 integerValue];
        }];
        NSString *key = [sortedKeys objectAtIndexedSubscript:row];
        LocusView *locusView = [self.stacksDocument.locusViews objectForKey:key];

        LocusCell* locusCell = (LocusCell *) cellView;
        locusCell.locusId.stringValue = locusView.locusId;
        locusCell.propertyField.stringValue = [NSString stringWithFormat:@"Parents %d Progeny %d \nSNPS %d"
                ,0,locusView.genotypes.count,locusView.snps.count];
//        locusCell.consensusField.stringValue = locusView.consensus;


        // START FANCY
        NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:locusView.consensus];

        [string beginEditing];
        NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blueColor], NSForegroundColorAttributeName,
                [NSColor grayColor], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier Bold" size:14.0], NSFontAttributeName,
                nil];
        for (SnpView* snpView in locusView.snps) {
            NSRange selectedRange = NSMakeRange(snpView.column, 1);
            [string setAttributes:attributes range:selectedRange];
        }
        [string endEditing];


        locusCell.consensusField.attributedStringValue = string;

        return locusCell;
    }
    else if ([[tableView identifier] isEqualToString:@"StacksTableView"]) {
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
                    case 0:
                        cellView.textField.stringValue = stacksView.reference.sequence;
                        break;
                    case 1:
                        cellView.textField.stringValue = stacksView.consensus.sequence;
                        break;
                    case 2:
                        cellView.textField.stringValue = stacksView.model.sequence;
                        break;
                    default:
                        NSString *sequenceString = [(StackEntry *) [stacksView.stackEntries objectAtIndex:row - 3] sequence];
//                        NSMutableAttributedString *string = [self decorateSnps:sequenceString snps:stacksView.snps];
                        NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:sequenceString];

                        [string beginEditing];
                        NSNumber *snpIndex;
                        NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                                [NSColor blueColor], NSForegroundColorAttributeName,
                                [NSColor grayColor], NSBackgroundColorAttributeName,
                                [NSFont fontWithName:@"Courier Bold" size:14.0], NSFontAttributeName,
                                nil];
                        for (snpIndex in stacksView.snps) {
                            NSRange selectedRange = NSMakeRange([snpIndex intValue], 1);
                            [string setAttributes:attributes range:selectedRange];
                        }
                        [string endEditing];


                        cellView.textField.attributedStringValue = string;
                        cellView.textField.font = [NSFont fontWithName:@"Courier" size:14];
                }
            }
            else {
                NSLog(@"not sure what that column is %@", tableColumn.identifier);
                cellView.textField.stringValue = @"";
            }
        }
        return cellView;
    }
    else {
        NSLog(@"could not find table %@", [tableView identifier]);
        return cellView;
    }


}

- (NSTableCellView *)handleGenotypesTable:(NSString *)column row:(NSInteger)row cell:(NSTableCellView *)cellView {
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
//            GenotypeEntry *genotypeEntry = [locusView.progeny objectAtIndex:progenyIndex];
//            GenotypeEntry *genotypeEntry = [locusView.genotypes valueForKey:[NSString stringWithFormat:@"%ld", progenyIndex]];
            NSArray *sortedKeys = [[locusView.genotypes allKeys] sortedArrayUsingComparator:^(NSString * obj1,NSString * obj2){
                return [obj1 compare:obj2];
            }];
            NSString *key = [sortedKeys objectAtIndex:progenyIndex-1];
            GenotypeEntry *genotypeEntry = [locusView.genotypes valueForKey:key];
            cellView.textField.stringValue = [NSString stringWithFormat:@"%@  %@", genotypeEntry.name, [genotypeEntry render]];
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



- (void)genotypeSelected:(id)tableView {
//    NSTableCellView *tableCellView = [tableView selectedCell];
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
    NSLog(@"loading genotypes %d",locusView.genotypes.count);

    if (index + 1 < locusView.genotypes.count) {
//            GenotypeEntry *entry = (GenotypeEntry *) [locusView.progeny objectAtIndex:index + 1];
        NSString *key = [[locusView.genotypes allKeys] objectAtIndex:index+1];
        GenotypeEntry *genotypeEntry = [locusView.genotypes valueForKey:key];
        self.selectedGenotype = genotypeEntry ;

        NSLog(@"entry ID: %d",genotypeEntry.sampleId);

        NSLog(@"loading %@ tag - %d",genotypeEntry.name, genotypeEntry.tagId);
        StacksView *stacksView = [self.stacksLoader loadStacksView:genotypeEntry.name atPath:@"/tmp/stacks_tut/" forTag:genotypeEntry.tagId];
        self.selectedStacks = stacksView;


//        [self.genotypesScrollView setHidden:TRUE];
//        [self.snpsScrollView setHidden:FALSE];
//        [self.stacksTableView setHidden:FALSE];
//        [self.stacksTableView reloadData];




    }
    else {
        NSLog(@"invalid selection");
        self.selectedStacks = nil ;
    }
    [self showTagsTable:self.selectedStacks];
}

@end


