//
//  DatumArrayController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/17/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "DatumArrayController.h"

@implementation DatumArrayController


- (void)awakeFromNib {
    [super awakeFromNib];

    [self setSortDescriptors:[NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"sampleId" ascending:YES selector:@selector(compare:)]]];
}

@end
