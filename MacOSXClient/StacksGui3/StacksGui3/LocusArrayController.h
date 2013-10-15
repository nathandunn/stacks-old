//
//  LocusArrayController.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/17/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface LocusArrayController : NSArrayController{
}

@property NSInteger minSnpValue;
@property NSInteger maxSnpValue;
- (IBAction) setMinSnpValue:(NSInteger) sender;
- (IBAction) setMaxSnpValue:(NSInteger) sender;
@end

