//
//  LocusMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "LocusMO.h"
#import "DatumMO.h"
#import "LocusAlleleMO.h"
#import "LocusSnpMO.h"


@implementation LocusMO

@dynamic consensus;
@dynamic length;
@dynamic locusId;
@dynamic marker;
@dynamic ratio;
@dynamic alleles;
@dynamic datums;
@dynamic snps;


- (NSAttributedString *)renderConsensus{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.consensus];
    [string beginEditing];
//        NSNumber *snpIndex;
    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blueColor], NSForegroundColorAttributeName,
            [NSColor grayColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
            nil];
    for (LocusSnpMO *snp in self.snps) {
        NSRange selectedRange = NSMakeRange([snp.column unsignedIntegerValue], 1);
        [string setAttributes:attributes range:selectedRange];
    }
    [string endEditing];
    return string ;
}

@end
