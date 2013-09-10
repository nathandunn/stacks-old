//
//  ModelStackEntryMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "ModelStackEntryMO.h"


@implementation ModelStackEntryMO

- (NSAttributedString *)renderEntryId{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:@""];
    return string ;
}

- (NSAttributedString *)renderRelationship{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.relationship];
    return string ;
}

- (NSAttributedString *)renderSequence {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.sequence];
    return string;
}

@end
