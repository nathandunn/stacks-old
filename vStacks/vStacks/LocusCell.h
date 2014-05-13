//
// Copyright (c) 2014 University of Oregon.  All rights reserved.
// Created by Nathan Dunn on 3/18/13.
//
//
//


#import <Foundation/Foundation.h>

//@class LocusView;


@interface LocusCell : NSTableCellView{
}

@property(assign) IBOutlet NSTextField *locusId;
@property(assign) IBOutlet NSTextField *propertyField;
@property(assign) IBOutlet NSTextField *consensusField;

@end