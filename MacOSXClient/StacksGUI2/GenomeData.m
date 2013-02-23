//
//  GenomeData.m
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import "GenomeData.h"

@implementation GenomeData

- (id)initWithGene:(NSString*)gene letters:(NSString*)letters{
    if ((self = [super init])) {
        self.gene= gene;
        self.letters = letters;
    }
    return self;
}


@end
