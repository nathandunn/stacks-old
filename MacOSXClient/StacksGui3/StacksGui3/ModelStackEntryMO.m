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
    NSDictionary *blockAttribute= [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blackColor], NSForegroundColorAttributeName,
            [NSColor whiteColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
            nil];
    [string setAttributes:blockAttribute range:NSMakeRange(0, self.sequence.length)];
    return string;
}

@end
