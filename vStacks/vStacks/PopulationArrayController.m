//
//  PopulationArrayController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/17/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "PopulationArrayController.h"

@implementation PopulationArrayController

- (void)awakeFromNib {
    [super awakeFromNib];
    [self setSortDescriptors:[NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"populationId" ascending:YES selector:@selector(compare:)]]];

}



- (NSArray *)arrangedObjects {
    NSMutableArray *defaultArrangedObjects = [NSMutableArray arrayWithArray:[super arrangedObjects]];

    NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:
            @"All", @"name", @"No Filter", @"value", nil];
    [defaultArrangedObjects insertObject:dictionary atIndex:0];

    return defaultArrangedObjects;
}


@end
