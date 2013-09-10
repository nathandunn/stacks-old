//
// Created by Nathan Dunn on 9/9/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>


@interface StacksApplicationController : NSDocumentController

@property(weak) IBOutlet NSPanel *progressPanel ;
@property(weak) IBOutlet NSProgressIndicator *loadProgress;
- (IBAction) importDocument:(id)sender;


@end