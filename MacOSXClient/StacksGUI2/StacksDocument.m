//
//  StacksDocument.m
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import "StacksDocument.h"
#import "LocusView.h"

@implementation StacksDocument

- (id)initWithLocusData:(LocusView *)locusData {
    if ((self = [super init])) {
        self.locusData = locusData;
    }
    return self ;
}


- (id)initWithMarker:(NSString*)marker consensusSequence:(NSString*)conensusSequence {
    if ((self = [super init])) {
        self.locusData = [[LocusView alloc] initWithId:marker consensus:conensusSequence];
    }
    return self ; 
}

@end
