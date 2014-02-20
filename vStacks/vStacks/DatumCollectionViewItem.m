//
// Created by Nathan Dunn on 5/26/13.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//
//
//


#import "DatumCollectionViewItem.h"
#import "DatumView.h"

@implementation DatumCollectionViewItem {

}
-(void)setSelected:(BOOL)flag
{
    [super setSelected:flag];
    [(DatumView *)[self view] setSelected:flag];
    [(DatumView *)[self view] setNeedsDisplay:YES];

    [[self view] invalidateIntrinsicContentSize];
    
}


@end