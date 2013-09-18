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
#import "SampleMO.h"


@implementation LocusMO

@dynamic consensus;
@dynamic length;
@dynamic locusId;
@dynamic marker;
@dynamic ratio;
@dynamic alleles;
@dynamic datums;
@dynamic snps;

- (NSInteger) countProgeny{
    NSInteger count =0 ;
    for(DatumMO *datumMO in self.datums){
        NSString* sampleName = datumMO.sample.name ;
        if([sampleName rangeOfString:@"male"].location==NSNotFound){
            ++count ;
        }
    }
    return count  ;
}

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
