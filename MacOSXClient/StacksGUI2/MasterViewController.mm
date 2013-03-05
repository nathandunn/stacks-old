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

@interface MasterViewController ()

@property (weak) IBOutlet NSTableView *filesTableView;
@property (weak) IBOutlet NSTextField *headerField;
@property (weak) IBOutlet NSTextField *fastaField;

@end

@implementation MasterViewController


// -------------------------------------------------------------------------------
//	awakeFromNib:
// -------------------------------------------------------------------------------
- (void)awakeFromNib
{
    [verticalSplitView setDelegate:self];	// we want a chance to affect the vertical split view coverage
}

// -------------------------------------------------------------------------------
//	splitView:effectiveRect:effectiveRect:forDrawnRect:ofDividerAtIndex
// -------------------------------------------------------------------------------
- (NSRect)splitView:(NSSplitView *)splitView effectiveRect:(NSRect)proposedEffectiveRect forDrawnRect:(NSRect)drawnRect ofDividerAtIndex:(NSInteger)dividerIndex
{
    NSRect effectiveRect = drawnRect;

    if (splitView == verticalSplitView)
    {
        // don't steal as much from the scroll bar as NSSplitView normally would
        effectiveRect.origin.x -= 2.0;
        effectiveRect.size.width += 6.0;

    }

    return effectiveRect;
}

// -------------------------------------------------------------------------------
//	splitView:additionalEffectiveRectOfDividerAtIndex:dividerIndex:
// -------------------------------------------------------------------------------
- (NSRect)splitView:(NSSplitView *)splitView additionalEffectiveRectOfDividerAtIndex:(NSInteger)dividerIndex
{
    // we have a divider handle next to one of the split views in the window
    if (splitView == verticalSplitView)
        return [dividerHandleView convertRect:[dividerHandleView bounds] toView:splitView];
    else
        return NSZeroRect;
}

// -------------------------------------------------------------------------------
//	constrainMinCoordinate:proposedCoordinate:index
// -------------------------------------------------------------------------------
- (CGFloat)splitView:(NSSplitView *)splitView constrainMinCoordinate:(CGFloat)proposedCoordinate ofSubviewAt:(NSInteger)index
{
    CGFloat constrainedCoordinate = proposedCoordinate;
    if (splitView == verticalSplitView)
    {
        // the primary vertical split view is asking for a constrained size
        constrainedCoordinate = proposedCoordinate + 120.0;
    }
    else if (splitView == horizontalSplitView)
    {
        // the horizontal split view between mailboxes and activity view
        constrainedCoordinate = proposedCoordinate + 200.0;
    }

    return constrainedCoordinate;
}

- (NSInteger)numberOfRowsInTableView:(NSTableView *)tableView {
    return [self.stacksDocuments count];
}

- (NSView *)tableView:(NSTableView *)tableView viewForTableColumn:(NSTableColumn *)tableColumn row:(NSInteger)row {

    // Get a new ViewCell
    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];

    // we want data for the row . . . .
    NSArray *sortedKeys = [[self.stacksDocuments allKeys] sortedArrayUsingComparator:(NSComparator)^(id obj1,id obj2){
        return [obj1 integerValue] - [obj2 integerValue];
    }];
    NSString *key = [sortedKeys objectAtIndexedSubscript:row];
    StacksDocument *stacksDoc = [self.stacksDocuments objectForKey:key];

    // Since this is a single-column table view, this would not be necessary.
    // But it's a good practice to do it in order by remember it when a table is multicolumn.
    if( [tableColumn.identifier isEqualToString:@"IdColumn"] )
    {
        cellView.textField.stringValue = stacksDoc.locusData.locusId;
        return cellView;
    }
    else
    if( [tableColumn.identifier isEqualToString:@"SnpColumn"] )
    {
        NSMutableArray *snps = stacksDoc.locusData.snps;
        if([snps count]>0){
            cellView.textField.stringValue = [NSString stringWithFormat:@"Yes [%dnuc]",[snps count]];
        }
        else{
            cellView.textField.stringValue = @"None";
        }
        return cellView;
    }
//    else
//    if( [tableColumn.identifier isEqualToString:@"ConsensusColumn"] )
//    {
//        cellView.textField.stringValue = stacksDoc.locusData.locusId;
//        return cellView;
//    }
    else
    if( [tableColumn.identifier isEqualToString:@"ParentsColumn"] )
    {
//        cellView.textField.stringValue = [NSString stringWithFormat:@"%d",[stacksDoc.locusData matchingParents]];
        cellView.textField.integerValue = [stacksDoc.locusData matchingParents];
//        NSLog(@"in the PARENTS column! value set: %@",stacksDoc.locusData.locusId);
        return cellView;
    }
    else
    if( [tableColumn.identifier isEqualToString:@"ProgenyColumn"] )
    {
        NSUInteger count = [[stacksDoc.locusData progeny] count];
        cellView.textField.stringValue = [NSString stringWithFormat:@"%d / %d",count,count];
        return cellView;
    }
    else
    if( [tableColumn.identifier isEqualToString:@"MarkerColumn"] )
    {
        cellView.textField.stringValue = stacksDoc.locusData.marker;
//        NSLog(@"MARKER column! value set: %@",stacksDoc.locusData.locusId);
        return cellView;
    }
    else
    if( [tableColumn.identifier isEqualToString:@"RatioColumn"] )
    {
        cellView.textField.stringValue = @"aa: 45 (51.7%) bb:42 (48.3%)";
        return cellView;
    }
    else
    if( [tableColumn.identifier isEqualToString:@"GenotypesColumn"] )
    {
        cellView.textField.integerValue = [stacksDoc.locusData genotypes] ;
        return cellView;
    }


    return cellView;
}

-(StacksDocument*) selectedDoc{
    NSInteger selectedRow = [self.filesTableView selectedRow];
    if(selectedRow >= 0 && self.stacksDocuments.count > selectedRow){
        NSString *key = [NSString stringWithFormat:@"%d",selectedRow];
        StacksDocument *stacksDocument = [self.stacksDocuments objectForKey:key];
        return stacksDocument;
    }
    return nil ;
}

-(void) setDetailInfo:(StacksDocument*) doc{
    NSString *name =@"";
    NSString *data =@"";

    if(doc!=nil){
        name = doc.locusData.locusId;
        data = doc.locusData.consensus;
    }

//    [self.testField setStringValue:test];
    [self.headerField setStringValue:name];
    [self.fastaField setStringValue:data];
}

- (void)tableViewSelectionDidChange:(NSNotification *)aNotification
{
    StacksDocument *selectedDoc = [self selectedDoc];

    // Update info
    [self setDetailInfo:selectedDoc];
}



@end



//- (id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil
//{
//    self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
//    if (self) {
//        // Initialization code here.
//    }
//
//    return self;
//}
//
//- (NSView *)tableView:(NSTableView *)tableView viewForTableColumn:(NSTableColumn *)tableColumn row:(NSInteger)row {
//
//    // Get a new ViewCell
//    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];
//
//
//    // Since this is a single-column table view, this would not be necessary.
//    // But it's a good practice to do it in order by remember it when a table is multicolumn.
//    if( [tableColumn.identifier isEqualToString:@"GeneColumn"] )
//    {
//        StacksDocument *bugDoc = [self.data objectAtIndex:row];
////        cellView.imageView.image = bugDoc.thumbImage;
//        cellView.textField.stringValue = bugDoc.locusData.locusId;
//        return cellView;
//    }
//    return cellView;
//}
//
//-(void) loadView{
//    [super loadView];
//    // customzie below this line
//}
//
//
//
//
//
//@end
