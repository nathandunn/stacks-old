//
//  PopulationMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "PopulationMO.h"
#import "SampleMO.h"


@implementation PopulationMO

@dynamic name;
@dynamic populationId;
@dynamic samples;
@dynamic metaData;

- (NSString *)annotatedName {

    NSNumberFormatter *numberFormatter = [[NSNumberFormatter alloc] init];
    numberFormatter.numberStyle = NSNumberFormatterNoStyle;

    NSNumber* number = [numberFormatter numberFromString:self.name];
    if(number==nil){
        return self.name ;
    }
    else{
        return [NSString stringWithFormat:@"Population %@",self.name];
    }
}
@end
