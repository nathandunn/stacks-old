//
//  ReferenceStackEntryMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "ReferenceStackEntryMO.h"


@implementation ReferenceStackEntryMO

- (NSAttributedString *)renderEntryId{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:@""];
    return string ;
}

- (NSAttributedString *)renderRelationship{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:@""];
    return string ;
}

- (NSAttributedString *)renderSequence {
    // create a string from 0-9 for sequenceSize
    NSUInteger sequenceSize = self.sequence.length ;

    NSString *string = [[NSString alloc] init];

    for (NSUInteger i = 0; i < sequenceSize; i++) {
        string = [string stringByAppendingFormat:@"%ld", i % 10];
    }

    NSMutableAttributedString *attributedString = [[NSMutableAttributedString alloc] initWithString:string];

    NSUInteger nextCount = 10;
    NSUInteger i = 0;
    NSUInteger next = 10;
    bool highlight1 = true;


    [attributedString beginEditing];

    while (i < sequenceSize) {
        next = i + nextCount;
        if (next > sequenceSize) {
            next = sequenceSize;
        }

//        NSRange range =
        NSRange selectedRange = NSMakeRange(i, next - i);

        NSDictionary *attributes;
        if (highlight1) {
            attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor blueColor], NSForegroundColorAttributeName,
                    [NSColor grayColor], NSBackgroundColorAttributeName,
                    [NSFont fontWithName:@"Courier" size:14], NSFontAttributeName,
                    nil];
            highlight1 = false;
        }
        else {
            attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor grayColor], NSForegroundColorAttributeName,
                    [NSColor blueColor], NSBackgroundColorAttributeName,
                    [NSFont fontWithName:@"Courier" size:14], NSFontAttributeName,
                    nil];
            highlight1 = true;

        }
        [attributedString setAttributes:attributes range:selectedRange];
        i += nextCount;

    }
    [attributedString endEditing];


    return attributedString;
}

@end
