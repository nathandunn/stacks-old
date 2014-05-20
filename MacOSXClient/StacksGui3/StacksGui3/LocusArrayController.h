//
//  LocusArrayController.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/17/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@class LocusRepository;

@interface LocusArrayController : NSArrayController{
}

@property NSInteger minSnpValue;
@property NSInteger maxSnpValue;
@property (nonatomic, retain) NSString*chromosomeLocation;
@property double minBasePairs;
@property double maxBasePairs;
@property(nonatomic, strong) LocusRepository* locusRepository;


- (IBAction) setMinSnpValue:(NSInteger) sender;
- (IBAction) setMaxSnpValue:(NSInteger) sender;

- (IBAction)writeLocationValue:(id)sender;
- (IBAction)writeMinBasePairs:(id)sender;
- (IBAction)writeMaxBasePairs:(id)sender;
@end
