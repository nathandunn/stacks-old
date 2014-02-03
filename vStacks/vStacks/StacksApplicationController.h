//
// Created by Nathan Dunn on 9/9/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class ProgressController ;

@interface StacksApplicationController : NSDocumentController{
//    ProgressController *progressController ;
}

@property(nonatomic, strong) NSNumberFormatter *numberFormatter;

- (IBAction) importDocument:(id)sender;


@end

