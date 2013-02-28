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
//@property (weak) IBOutlet NSTextField *testField;
@property (weak) IBOutlet NSTextField *headerField;
@property (weak) IBOutlet NSTextField *fastaField;
//@property (weak) IBOutlet NSTextField *markerField;


//@property (weak) IBOutlet NSTextField *testField ;
//@property (weak) IBOutlet NSTextField *headerField ;
//@property (weak) IBOutlet NSTextField *fastaField ;

@end

@implementation MasterViewController

- (id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil
{
    self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
    if (self) {
        // Initialization code here.
    }
    
    return self;
}

- (NSView *)tableView:(NSTableView *)tableView viewForTableColumn:(NSTableColumn *)tableColumn row:(NSInteger)row {
    
    // Get a new ViewCell
    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];
    
    
    // Since this is a single-column table view, this would not be necessary.
    // But it's a good practice to do it in order by remember it when a table is multicolumn.
    if( [tableColumn.identifier isEqualToString:@"GeneColumn"] )
    {
        StacksDocument *bugDoc = [self.data objectAtIndex:row];
//        cellView.imageView.image = bugDoc.thumbImage;
        cellView.textField.stringValue = bugDoc.locusData.locusId;
        return cellView;
    }
    return cellView;
}

-(void) loadView{
    [super loadView];
    // customzie below this line
}


-(StacksDocument*) selectedDoc{
    NSInteger selectedRow = [self.filesTableView selectedRow];
    if(selectedRow >= 0 && self.data.count > selectedRow){
        StacksDocument *stacksDocument = [self.data objectAtIndex:selectedRow];
        return stacksDocument;
    }
    return nil ;
}

-(void) setDetailInfo:(StacksDocument*) doc{
    NSString *name =@"";
    NSString *data =@"";
    NSString *test =@"";
    
    if(doc!=nil){
        name = doc.locusData.locusId;
        data = doc.locusData.consensus;
        test = @"testis working I think";
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



- (NSInteger)numberOfRowsInTableView:(NSTableView *)tableView {
    return [self.data count];
}

@end
