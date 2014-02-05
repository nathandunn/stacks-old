//
//  DatumView.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/26/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@class DatumMO;

//@interface DatumView : NSCollectionViewItem
@interface DatumView : NSView{
    BOOL selected ;
}

@property (readwrite) BOOL selected;
@property (weak) IBOutlet NSTextField* haplotypeField;



@end
