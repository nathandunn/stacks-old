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

- (id)initWithLocusView:(LocusView *)locusData {
    if ((self = [super init])) {
        self.locusData = locusData;
        self.locusId = locusData.locusId;
    }
    return self ;
}


@end

