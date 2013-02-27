//
//  StacksDocument.m
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import "StacksDocument.h"
#import "GenomeData.h"

@implementation StacksDocument


- (id)initWithMarker:(NSString*)marker consensusSequence:(NSString*)conensusSequence {
    if ((self = [super init])) {
        self.data = [[GenomeData alloc] initWithMarker: marker consensusSequence:conensusSequence];
    }
    return self ; 
}

@end
