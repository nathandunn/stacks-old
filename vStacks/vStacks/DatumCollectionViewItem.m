//
// Created by Nathan Dunn on 5/26/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "DatumCollectionViewItem.h"
#import "DatumView.h"

@implementation DatumCollectionViewItem {

}
-(void)setSelected:(BOOL)flag
{
    [super setSelected:flag];
    [(DatumView *)[self view] setSelected:flag];
//    [(DatumView *)[self view] setNeedsDisplay:YES];

    [[self view] invalidateIntrinsicContentSize];
    
}


@end