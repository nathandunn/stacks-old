//
//  StacksDocument.h
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>

@class LocusView;

@interface StacksDocument : NSObject

@property (strong) LocusView *locusData;

// rename to init with locaus data
- (id)initWithLocusData:(LocusView *)locusData;
- (id)initWithMarker:(NSString*)marker consensusSequence:(NSString*)consensusSequence ;


@end
