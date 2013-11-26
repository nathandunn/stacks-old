//
//  ConsensusStackEntryMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "ConsensusStackEntryMO.h"
#import "DatumSnpMO.h"
#import "DatumMO.h"
#import "LocusMO.h"
#import "LocusSnpMO.h"


@implementation ConsensusStackEntryMO

- (NSAttributedString *)renderEntryId {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:@""];
    return string;
}

- (NSAttributedString *)renderRelationship {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.relationship];
    return string;
}

- (NSAttributedString *)renderSequence {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.sequence];
    [string beginEditing];

    NSDictionary *blockAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blackColor], NSForegroundColorAttributeName,
//            [NSColor whiteColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
            nil];
    [string setAttributes:blockAttribute range:NSMakeRange(0, self.sequence.length)];

    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor redColor], NSForegroundColorAttributeName,
            [NSColor grayColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier Bold" size:14.0], NSFontAttributeName,
            nil];


//    NSSet *locusSnps = self.datum.locus.snps ;
    // process locus snps
    NSMutableArray *locusSnpColumns = [self getLocusSnps];

    for (NSNumber *snpColumn in locusSnpColumns) {
        NSRange selectedRange = NSMakeRange([snpColumn unsignedIntegerValue], 1);
        [string setAttributes:attributes range:selectedRange];
    }
    [string endEditing];
    return string;
}

@end
