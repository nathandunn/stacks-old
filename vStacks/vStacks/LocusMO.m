//
//  LocusMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "LocusMO.h"
#import "DatumMO.h"
#import "SampleMO.h"
#import "ColorGenerator.h"


@implementation LocusMO

@dynamic basePairs;
@dynamic parentCount;
@dynamic progenyCount;
@dynamic chromosome;
@dynamic type;
@dynamic strand;
@dynamic consensus;
@dynamic length;
@dynamic locusId;
@dynamic marker;
@dynamic ratio;
@dynamic alleleData;
//@dynamic datums;
@dynamic snpData;
@dynamic metaData;

@synthesize colorGenerator;

- (id)init {
    self = [super init];
    if (self) {
    }

    return self;
}


- (NSInteger)countParents {
    NSInteger count = 0;
    // TODO: fix

    //    for (DatumMO *datumMO in self.datums) {
//        NSString *sampleName = datumMO.sample.name;
//        if ([sampleName rangeOfString:@"male"].location != NSNotFound) {
//            ++count;
//        }
//    }
    return count;
}


- (NSUInteger)countProgeny {
    NSUInteger count = 0;
    // TODO: fix
//    NSLog(@"number of datums!! %ld",self.datums.count);
//    for (DatumMO *datumMO in self.datums) {
//        NSString *sampleName = datumMO.sample.name;
////        NSLog(@"sample name %@",sampleName);
//        if ([sampleName rangeOfString:@"male"].location == NSNotFound) {
////            NSLog(@"not found!!") ;
//            ++count;
//        }
//    }
    return count;
}

- (NSAttributedString *)renderChromosome {
    NSString *inputString = @"";
//    NSLog(@"self.chromosome %@", self.chromosome);
    if (self.chromosome != nil && self.chromosome.length > 0) {
        inputString = [NSString stringWithFormat:@"%@ %@ Mb %@", self.chromosome, self.basePairs, self.strand];
    }
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:inputString];
    return string;
}

- (NSAttributedString *)renderDescription {
    NSNumber *parentCount = self.parentCount;
    NSUInteger progenyCount = self.progenyCount.unsignedIntegerValue ;
    // TODO: convert

    NSError *error;
    NSDictionary *json = [NSJSONSerialization JSONObjectWithData:self.snpData options:kNilOptions error:&error];
    NSUInteger snpCount = json.count;

    NSString *chromosomeString = @"";
    if (self.chromosome != nil && self.chromosome.length > 0) {
        double basePairsValue = [self.basePairs doubleValue] / 1000000.0f;
        NSString *basePairsString = [NSString stringWithFormat:@"%1.2f", basePairsValue];
        chromosomeString = [NSString stringWithFormat:@"%@ %@Mb %@", self.chromosome, basePairsString, self.strand];
    }

    NSString *inputString;
//    NSLog(@"type %@",self.type);
    if ([self.type isEqualToString:@"GeneticMap"]) {
//        parentString = @"Parents";
        inputString = [NSString stringWithFormat:@"Parents %ld Prog %ld Snps %ld\n%@", parentCount.integerValue, progenyCount-parentCount.unsignedIntegerValue, snpCount, chromosomeString];
    }
    else {
//        parentString = @"Samples";
        inputString = [NSString stringWithFormat:@"Samples %ld Snps %ld\n%@", progenyCount, snpCount, chromosomeString];
    }


    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:inputString];
    [string beginEditing];

    NSMutableParagraphStyle *paragraph = [[NSMutableParagraphStyle alloc] init];
    paragraph.alignment = NSRightTextAlignment;
    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blackColor], NSForegroundColorAttributeName,
            [NSFont fontWithName:@"Helvetica" size:10.0], NSFontAttributeName,
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

    // TODO: convert
    NSError *error;
    NSDictionary *snpJson = [NSJSONSerialization JSONObjectWithData:self.snpData options:kNilOptions error:&error];
    for (NSDictionary *snp in snpJson) {
        NSInteger startRange = [[snp valueForKey:@"column"] integerValue];
        NSRange selectedRange = NSMakeRange(startRange, 1);
        [string addAttributes:attributes range:selectedRange];
    }


    [string endEditing];
    return string;
}


- (ColorGenerator *)getColorGenerator {
    if (colorGenerator == nil) {
        colorGenerator = [[ColorGenerator alloc] init];
    }
    return colorGenerator;
}
@end
