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

- (id)initWithGene:(NSString*)gene letters:(NSString*)letters {
    if ((self = [super init])) {
        self.data = [[GenomeData alloc] initWithGene:gene letters:letters];
    }
    return self ; 
}

@end
