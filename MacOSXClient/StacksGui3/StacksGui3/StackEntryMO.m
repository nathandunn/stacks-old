//
//  StackEntryMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StackEntryMO.h"
#import "DatumMO.h"
#import "SnpMO.h"
#import "DatumSnpMO.h"
#import "LocusSnpMO.h"
#import "LocusMO.h"
#import "ConsensusStackEntryMO.h"
#import "ColorGenerator.h"


@implementation StackEntryMO

@dynamic block;
@dynamic entryId;
@dynamic relationship;
@dynamic sequence;
@dynamic sequenceId;
@dynamic datum;

@synthesize colorGenerator;


- (NSAttributedString *)renderEntryId {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:[NSString stringWithFormat:@"%@", self.entryId]];
    return string;
}

- (NSAttributedString *)renderRelationship {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:[NSString stringWithFormat:@"%@", self.relationship]];
    [string beginEditing];
    NSDictionary *attributes;

    if ([self.relationship isEqualToString:@"primary"]) {
        attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor greenColor], NSForegroundColorAttributeName,
                nil];
    }
            // secondary
    else {
        attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor redColor], NSForegroundColorAttributeName,
                nil];
    }
    [string setAttributes:attributes range:NSMakeRange(0, self.relationship.length)];
    [string endEditing];
    return string;
}

- (NSAttributedString *)renderSequence {
//    NSMutableAttributedString *attributedString = [[NSMutableAttributedString alloc] initWithString:self.sequence];

    NSMutableAttributedString *sequenceAttributedString = [[NSMutableAttributedString alloc] initWithString:self.sequence];

    [sequenceAttributedString beginEditing];
    NSDictionary *blockAttribute;

    // for both a catalog and local snp
    NSDictionary *snpAttribute;

    // for both a catalog and local snp
    NSDictionary *catalogSnpOnlyAttribute;

    NSDictionary *defectAttribute;


    /**
     https://casspr.fogbugz.com/default.asp?1550#11745
     Catalog SNP but not a local SNP, column (background) color is blueish: #d8f3ff
     Catalog SNP and local SNP, column (background) color is greyish: #ccddd4
     Catalog SNP and local SNP, column (foreground) color is #a93535, bold for first allele, #29356c, bold for second allele
     
     SNP background color is #EEE
     */



    // DEFINE BLOCK
    NSInteger blockValue = [self.block integerValue];

    NSLog(@"START defining attributes!! %@", self.getColorGenerator);


    defectAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blackColor], NSForegroundColorAttributeName,
            [self.getColorGenerator colorWithHexString:@"fcdba7"], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
            nil];

    if (blockValue % 2 == 0) {
//        NSLog(@"block 1->%@",self.block);
        catalogSnpOnlyAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blackColor], NSForegroundColorAttributeName,
                [self.getColorGenerator colorWithHexString:@"d8f3ff"], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier-Bold" size:14.0], NSFontAttributeName,
                nil];

        snpAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [self.getColorGenerator colorWithHexString:@"a93535"], NSForegroundColorAttributeName,
//                [self.getColorGenerator colorWithHexString:@"ccddd4"], NSBackgroundColorAttributeName,
                [self.getColorGenerator colorWithHexString:@"d8f3ff"], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier-Bold" size:14.0], NSFontAttributeName,
                nil];
        blockAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blackColor], NSForegroundColorAttributeName,
                [NSColor whiteColor], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
                nil];
    }
    else {
//        NSLog(@"block 0->%@",self.block);

        catalogSnpOnlyAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blackColor], NSForegroundColorAttributeName,
                [self.getColorGenerator colorWithHexString:@"d8f3ff"], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier-Bold" size:14.0], NSFontAttributeName,
                nil];

        snpAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [self.getColorGenerator colorWithHexString:@"29356c"], NSForegroundColorAttributeName,
//                [self.getColorGenerator colorWithHexString:@"ccddd4"], NSBackgroundColorAttributeName,
                [self.getColorGenerator colorWithHexString:@"d8f3ff"], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier-Bold" size:14.0], NSFontAttributeName,
                nil];
        blockAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blackColor], NSForegroundColorAttributeName,
                [self.getColorGenerator colorWithHexString:@"dddddd"], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
                nil];
    }

    NSLog(@"END defining attributes!!");

    [sequenceAttributedString setAttributes:blockAttribute range:NSMakeRange(0, self.sequence.length)];

    // get a local snp
    NSSet *datumSnps = self.datum.snps;
    NSSet *locusSnps = self.datum.locus.snps;

    NSMutableArray *locusSnpColumns = [[NSMutableArray alloc] init];
    NSMutableArray *datumSnpColumns = [[NSMutableArray alloc] init];

    // process locus snps
    for (LocusSnpMO *locusSnp in locusSnps) {

        [locusSnpColumns addObject:locusSnp.column];
    }

    for (LocusSnpMO *datumSnp in datumSnps) {
        [datumSnpColumns addObject:datumSnp.column];
    }

    // color Snps
    for (NSNumber *column in locusSnpColumns) {
        NSRange selectedRange = NSMakeRange([column unsignedIntegerValue], 1);

        if ([datumSnpColumns containsObject:column]) {
            [sequenceAttributedString setAttributes:snpAttribute range:selectedRange];
        }
        else {
            [sequenceAttributedString setAttributes:catalogSnpOnlyAttribute range:selectedRange];
        }
    }


    // color any part of the consensus sequence that does not match.
    NSString *consensusSequence = self.datum.consensus.sequence;
    NSString *mySequence = [sequenceAttributedString string];
    if (![consensusSequence isEqualToString:mySequence]) {
        for (int i = 0; i < mySequence.length && i < consensusSequence.length; i++) {
            if ([consensusSequence characterAtIndex:i] != [mySequence characterAtIndex:i]) {
                // if is actually a defined SNP, then we ignore
                if (![locusSnpColumns containsObject:[NSNumber numberWithInt:i]]) {
                    NSRange selectedRange = NSMakeRange(i, 1);
                    [sequenceAttributedString setAttributes:defectAttribute range:selectedRange];
                }
            }
        }
    }

    [sequenceAttributedString endEditing];
    return sequenceAttributedString;

}

- (ColorGenerator *)getColorGenerator {
    if (colorGenerator == nil) {
        colorGenerator = [[ColorGenerator alloc] init];
    }
    return colorGenerator;
}
@end
