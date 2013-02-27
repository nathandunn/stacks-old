//
//  GenomeData.m
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import "GenomeData.h"

@implementation GenomeData

- (id)initWithMarker:(NSString *)marker consensusSequence:(NSString *)consensusSequence{
    if ((self = [super init])) {
        self.marker= marker;
        self.consensusSequence = consensusSequence;
    }
    return self;
}


@end
