//
//  GenotypesViewController.m
//  StacksGUI2
//
//  Created by Nathan Dunn on 3/5/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "GenotypesViewController.h"

@interface GenotypesViewController ()

@property (weak) IBOutlet NSTableView *genotypesTable;
@end

@implementation GenotypesViewController

- (id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil
{
    self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
    if (self) {
        // Initialization code here.
        _genotypesTable.dataSource = self;
        _genotypesTable.delegate = self;
    }
    
    return self;
}

- (NSInteger)numberOfRowsInTableView:(NSTableView *)tableView {
//    return 0;
    return 5;
}

- (NSView *)tableView:(NSTableView *)tableView viewForTableColumn:(NSTableColumn *)tableColumn row:(NSInteger)row {
//    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];
    NSTableCellView *cellView = [tableView makeViewWithIdentifier:tableColumn.identifier owner:self];
    return cellView;
}





@end
