//
//  DatumView.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/26/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "DatumView.h"

@implementation DatumView

@synthesize selected;



- (id)initWithFrame:(NSRect)frame
{
    self = [super initWithFrame:frame];
    if (self) {
        // Initialization code here.
    }

    return self;
}

- (void)drawRect:(NSRect)dirtyRect
{
    // Drawing code here.
    if (selected) {
//        [[NSColor blueColor] set];
        [[NSColor lightGrayColor] set];
        NSRectFill([self bounds]);
    }
}

@end
