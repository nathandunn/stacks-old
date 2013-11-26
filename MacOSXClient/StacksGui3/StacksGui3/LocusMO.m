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
#import "DepthMO.h"
#import "HaplotypeMO.h"
#import "ColorGenerator.h"


@implementation LocusMO

@dynamic basePairs;
@dynamic parentCount;
@dynamic chromosome;
@dynamic type;
@dynamic strand;
@dynamic consensus;
@dynamic length;
@dynamic locusId;
@dynamic marker;
@dynamic ratio;
@dynamic alleles;
@dynamic datums;
@dynamic snps;

@synthesize haplotypeOrder;

@synthesize colorGenerator;

- (id)init {
    self = [super init];
    if (self) {
    }

    return self;
}


- (NSInteger)countParents {
    NSInteger count = 0;
    for (DatumMO *datumMO in self.datums) {
        NSString *sampleName = datumMO.sample.name;
        if ([sampleName rangeOfString:@"male"].location != NSNotFound) {
            ++count;
        }
    }
    return count;
}

- (NSUInteger)countProgeny {
    NSUInteger count = 0;
//    NSLog(@"number of datums!! %ld",self.datums.count);
    for (DatumMO *datumMO in self.datums) {
        NSString *sampleName = datumMO.sample.name;
//        NSLog(@"sample name %@",sampleName);
        if ([sampleName rangeOfString:@"male"].location == NSNotFound) {
//            NSLog(@"not found!!") ;
            ++count;
        }
    }
    return count;
}

- (NSAttributedString *)renderChromosome {
    NSString *inputString = @"";
    NSLog(@"self.chromosome %@", self.chromosome);
    if (self.chromosome != nil && self.chromosome.length > 0) {
        inputString = [NSString stringWithFormat:@"%@ %@ Mb %@", self.chromosome, self.basePairs, self.strand];
    }
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:inputString];
    return string;
}

- (NSAttributedString *)renderDescription {
    NSNumber *parentCount = self.parentCount;
    NSUInteger progenyCount = self.countProgeny;
    NSUInteger snpCount = self.snps.count;
//    NSString *parentString;

    NSString *chromosomeString = @"";
    if (self.chromosome != nil && self.chromosome.length > 0) {
        double basePairsValue = [self.basePairs doubleValue] / 1000000.0f;
        NSString *basePairsString = [NSString stringWithFormat:@"%1.2f", basePairsValue];
        chromosomeString = [NSString stringWithFormat:@"%@ %@ Mb %@", self.chromosome, basePairsString, self.strand];
    }

    NSString *inputString ;
    NSLog(@"type %@",self.type);
    if ([self.type isEqualToString:@"GeneticMap"]) {
//        parentString = @"Parents";
        inputString = [NSString stringWithFormat:@"Parents %ld Prog %ld Snps %ld %@",  parentCount.integerValue, progenyCount, snpCount, chromosomeString];
    }
    else {
//        parentString = @"Samples";
        inputString = [NSString stringWithFormat:@"Samples %ld Snps %ld %@", progenyCount, snpCount, chromosomeString];
    }




    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:inputString];
    [string beginEditing];

    NSMutableParagraphStyle *paragraph = [[NSMutableParagraphStyle alloc] init];
    paragraph.alignment = NSRightTextAlignment;
    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blackColor], NSForegroundColorAttributeName,
//            [NSColor whiteColor], NSBackgroundColorAttributeName,
//            [NSFont fontWithName:@"Courier" size:12.0], NSFontAttributeName,
                    [NSFont fontWithName:@"Helvetica" size:12.0], NSFontAttributeName,
            paragraph, NSParagraphStyleAttributeName,
            nil];


    // Add attribute NSParagraphStyleAttributeName
//    [string addAttribute:NSParagraphStyleAttributeName
//                             value:paragraph range:NSMakeRange(0, string.length)];

    NSRange selectedRange = NSMakeRange(0, string.length);
    [string setAttributes:attributes range:selectedRange];
    [string endEditing];
    return string;
}

- (NSAttributedString *)renderConsensus {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.consensus];
    [string beginEditing];
//        NSNumber *snpIndex;

    NSRange range = NSMakeRange(0, string.length);
    float fontSize = 12;
    NSDictionary *globalAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blackColor], NSForegroundColorAttributeName,
//            [NSColor whiteColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier" size:fontSize], NSFontAttributeName,
            nil];
    [string setAttributes:globalAttribute range:range];


    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [self.getColorGenerator colorWithHexString:@"a93535"], NSForegroundColorAttributeName,
            [self.getColorGenerator colorWithHexString:@"ccddd4"], NSBackgroundColorAttributeName,
//                    [NSColor redColor], NSForegroundColorAttributeName,
//                    [NSColor grayColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier-Bold" size:fontSize], NSFontAttributeName,
            nil];
    for (LocusSnpMO *snp in self.snps) {
        NSRange selectedRange = NSMakeRange([snp.column unsignedIntegerValue], 1);
        [string addAttributes:attributes range:selectedRange];
    }


    [string endEditing];
    return string;
}

- (NSUInteger)lookupHaplotypeOrder:(NSString *)haplotype {
    if (haplotypeOrder == nil) {
        haplotypeOrder = [[NSMutableDictionary alloc] init];
        [haplotypeOrder setValue:[NSNumber numberWithInt:0] forKey:haplotype];
        return 0;
    }

    NSNumber *returnType = [haplotypeOrder objectForKey:haplotype];
    if (returnType == nil) {
        NSUInteger maxValue = 0;
        for (NSNumber *aValue in haplotypeOrder.allValues) {
            if ([aValue unsignedIntegerValue] > maxValue) {
                maxValue = [aValue unsignedIntegerValue];
            }
        }
        ++maxValue;
        [haplotypeOrder setValue:[NSNumber numberWithInt:maxValue] forKey:haplotype];
        // actually we get he max value here first .

        return maxValue;
    }
    else {
        return [returnType unsignedIntegerValue];
    }

}

- (NSUInteger)lookupDepthOrder:(DepthMO *)depthMO {
    // assume that it is there . . .
    NSNumber *order = depthMO.order;
    for (HaplotypeMO *haplotypeMO in depthMO.datum.haplotypes.allObjects) {
        if ([haplotypeMO.order isEqualToNumber:order]) {
            return [self lookupHaplotypeOrder:haplotypeMO.haplotype];
        }
    }

    return 0;
}

- (ColorGenerator *)getColorGenerator {
    if (colorGenerator == nil) {
        colorGenerator = [[ColorGenerator alloc] init];
    }
    return colorGenerator;
}
@end
