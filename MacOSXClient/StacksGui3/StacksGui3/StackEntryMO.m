//
//  StackEntryMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StackEntryMO.h"
#import "DatumMO.h"


@implementation StackEntryMO

@dynamic block;
@dynamic entryId;
@dynamic relationship;
@dynamic sequence;
@dynamic sequenceId;
@dynamic datum;


- (NSAttributedString *)renderEntryId{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:[NSString stringWithFormat:@"%@",self.entryId] ];
    return string ;
}

- (NSAttributedString *)renderRelationship{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:[NSString stringWithFormat:@"%@",self.relationship] ];
    [string beginEditing];
    NSDictionary *attributes ;

    if([self.relationship isEqualToString:@"primary"]){
    attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor greenColor], NSForegroundColorAttributeName,
            nil];
    }
    // secondary
    else{
        attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor redColor], NSForegroundColorAttributeName,
                nil];
    }
    [string setAttributes:attributes range:NSMakeRange(0, self.relationship.length)];
    [string endEditing];
    return string ;
}

- (NSAttributedString *)renderSequence {
//    NSMutableAttributedString *attributedString = [[NSMutableAttributedString alloc] initWithString:self.sequence];

    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.sequence];

    [string beginEditing];
//    NSNumber *snpIndex;
    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blueColor], NSForegroundColorAttributeName,
            [NSColor grayColor], NSBackgroundColorAttributeName,
//            [NSFont fontWithName:@"Courier Bold" size:14.0], NSFontAttributeName,
            [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
            nil];
//    for (snpIndex in self.datum.snps) {
//        NSRange selectedRange = NSMakeRange([snpIndex intValue], 1);
//        [string setAttributes:attributes range:selectedRange];
//    }
    [string setAttributes:attributes range:NSMakeRange(0, self.sequence.length)];
    [string endEditing];
    return string;


//    return attributedString ;
}
@end
