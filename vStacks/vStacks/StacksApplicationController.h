//
// Created by Nathan Dunn on 9/9/13.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//
//
//


#import <Foundation/Foundation.h>
n @class ProgressController ;

@interface StacksApplicationController : NSDocumentController{
//    ProgressController *progressController ;
}

@property(nonatomic, strong) NSNumberFormatter *numberFormatter;

- (IBAction) importDocument:(id)sender;
- (IBAction)provideFeedback:(id)sender;
- (IBAction) helpStacks:(id)sender;
- (IBAction) helpVStacks:(id)sender;
- (IBAction) license:(id)sender;


@end

