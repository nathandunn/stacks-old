//
//  PopulationArrayController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/17/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "PopulationArrayController.h"

@implementation PopulationArrayController

- (void)awakeFromNib {
    [super awakeFromNib];
    [self setSortDescriptors:[NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"populationId" ascending:YES selector:@selector(compare:)]]];

}

@end
