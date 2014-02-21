//
//  DatumMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "DatumMO.h"
//#import "ConsensusStackEntryMO.h"
//#import "DatumAlleleMO.h"
//#import "DatumSnpMO.h"
//#import "DepthMO.h"
//#import "HaplotypeMO.h"
//#import "LocusMO.h"
//#import "ModelStackEntryMO.h"
//#import "ReferenceStackEntryMO.h"
//#import "SampleMO.h"
//#import "StackEntryMO.h"
#import "ColorGenerator.h"


@implementation DatumMO

@dynamic name;
@dynamic sampleId;
@dynamic tagId;
//@dynamic alleles;
@dynamic alleleData;
//@dynamic consensus;
//@dynamic depths;
@dynamic depthData;
//@dynamic haplotypes;
@dynamic haplotypeData;
@dynamic locus;
//@dynamic model;
//@dynamic reference;
@dynamic sample;
//@dynamic snps;
@dynamic snpData;
//@dynamic stackEntries;
@dynamic stackData;
@dynamic metaData;

@synthesize colorGenerator;


- (NSAttributedString *)renderName {
    NSString *formattedString = [self.name stringByReplacingOccurrencesOfString:@"_" withString:@" "];
    formattedString = [formattedString capitalizedString];

    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:formattedString];

    NSMutableParagraphStyle *mutParaStyle = [[NSMutableParagraphStyle alloc] init];
    [mutParaStyle setAlignment:NSCenterTextAlignment];
    NSRange selectedRange = NSMakeRange(0, [[string string] length]);
    [string addAttributes:[NSDictionary dictionaryWithObject:mutParaStyle
                                                      forKey:NSParagraphStyleAttributeName] range:selectedRange];
//    NSDictionary *fontAttributes = [NSDictionary dictionaryWithObjectsAndKeys:[NSFont fontWithName:@"Courier" size:14], NSFontAttributeName, nil];
//    [string addAttributes:fontAttributes range:selectedRange];

    return string;
}

- (NSAttributedString *)renderHaplotypes {

    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] init];


    int i = 0;

    // TODO: convert
    NSError *error;
    NSDictionary *haplotypeJson = [NSJSONSerialization JSONObjectWithData:self.haplotypeData options:kNilOptions error:&error];
    NSUInteger haplotypeCount = haplotypeJson.count;
    for (NSDictionary *haplotype in haplotypeJson) {
        NSUInteger order = [[haplotype valueForKey:@"order"] unsignedIntegerValue];
        NSString *haplotypeString = [haplotype valueForKey:@"haplotype"];


        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:haplotypeString];
        NSRange selectedRange = NSMakeRange(0, haplotypeString.length);
        [appendString beginEditing];
        NSDictionary *attributes = [self generateColorForOrder:order];
        [appendString setAttributes:attributes range:selectedRange];
        NSDictionary *fontAttributes = [NSDictionary dictionaryWithObjectsAndKeys:[NSFont fontWithName:@"Courier" size:14], NSFontAttributeName, nil];
        [appendString addAttributes:fontAttributes range:selectedRange];
        [appendString endEditing];

        [string appendAttributedString:appendString];
//        if (i < self.haplotypes.count - 1) {
        if (i < haplotypeCount - 1) {
            [string appendAttributedString:[[NSAttributedString alloc] initWithString:@" / "]];
        }
        ++i;
    }

    NSMutableParagraphStyle *mutParaStyle = [[NSMutableParagraphStyle alloc] init];
    [mutParaStyle setAlignment:NSCenterTextAlignment];
    [string addAttributes:[NSDictionary dictionaryWithObject:mutParaStyle
                                                      forKey:NSParagraphStyleAttributeName]
                    range:NSMakeRange(0, [[string string] length])];

    return string;
}

- (NSAttributedString *)renderDepths {

    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] init];

    int i = 0;
    NSError *error;
    NSDictionary *haplotypeJson = [NSJSONSerialization JSONObjectWithData:self.haplotypeData options:kNilOptions error:&error];
    NSUInteger haplotypeCount = haplotypeJson.count;
    for (NSDictionary *haplotype in haplotypeJson) {
        NSUInteger order = [[haplotype valueForKey:@"order"] unsignedIntegerValue];
        NSString *depthString = [NSString stringWithFormat:@"%@", [haplotype valueForKey:@"depth"]];


        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:depthString];
        NSRange selectedRange = NSMakeRange(0, depthString.length);
        [appendString beginEditing];
        NSDictionary *attributes = [self generateColorForOrder:order];
        [appendString setAttributes:attributes range:selectedRange];
        NSDictionary *fontAttributes = [NSDictionary dictionaryWithObjectsAndKeys:[NSFont fontWithName:@"Courier" size:14], NSFontAttributeName, nil];
        [appendString addAttributes:fontAttributes range:selectedRange];
        [appendString endEditing];

        [string appendAttributedString:appendString];
        if (i < haplotypeCount - 1) {
            [string appendAttributedString:[[NSAttributedString alloc] initWithString:@" / "]];
        }
        ++i;
    }


    NSMutableParagraphStyle *mutParaStyle = [[NSMutableParagraphStyle alloc] init];
    [mutParaStyle setAlignment:NSCenterTextAlignment];
    [string addAttributes:[NSDictionary dictionaryWithObject:mutParaStyle
                                                      forKey:NSParagraphStyleAttributeName]
                    range:NSMakeRange(0, [[string string] length])];

    return string;
}


- (NSString *)generateColorStringForOrder:(NSUInteger)order {
    if (colorGenerator == nil) {
        colorGenerator = [[ColorGenerator alloc] init];
    }

    return [colorGenerator generateColorStringForOrder:order];
}


- (NSDictionary *)generateColorForOrder:(NSUInteger)order {
    if (colorGenerator == nil) {
        colorGenerator = [[ColorGenerator alloc] init];
    }

    NSColor *color = [colorGenerator generateColorForOrder:order];
    return [NSDictionary dictionaryWithObjectsAndKeys:
            color, NSForegroundColorAttributeName,
            nil];

}

- (NSMutableString *)renderHaplotypeHtml:(NSDictionary *)hapReorder{
    NSMutableString *returnHTML = [NSMutableString string];


    int i = 0;
    NSError *error;
    NSDictionary *haplotypeJson = [NSJSONSerialization JSONObjectWithData:self.haplotypeData options:kNilOptions error:&error];
    NSUInteger haplotypeCount = haplotypeJson.count;
    [returnHTML appendString:@"<div class='datum-depth'>"];
    for (NSDictionary *haplotype in haplotypeJson) {
//        NSUInteger order = [[haplotype valueForKey:@"order"] unsignedIntegerValue];
        NSNumber *order = [hapReorder objectForKey:[haplotype valueForKey:@"order"]];
//        NSLog(@"rednering order: %ld", order);
        NSString *haplotypeString = [NSString stringWithFormat:@"%@", [haplotype valueForKey:@"haplotype"]];

        NSMutableString *appendString = [NSMutableString stringWithString:haplotypeString];
        NSString *colorString = [self generateColorStringForOrder:order.unsignedIntegerValue];
        [appendString insertString:[NSString stringWithFormat:@"<font color='#%@'>", colorString] atIndex:0];
        [appendString appendString:@"</font>"];
//        NSDictionary *fontAttributes = [NSDictionary dictionaryWithObjectsAndKeys:[NSFont fontWithName:@"Courier" size:14], NSFontAttributeName, nil];

        [returnHTML appendString:appendString];
        if (i < haplotypeCount - 1) {
            [returnHTML appendString:@" / "];
        }
        ++i;
    }
    [returnHTML appendString:@"</div>"];

    return returnHTML;
}

- (NSMutableString *)renderDepthHtml:(NSDictionary *)hapReorder{
    NSMutableString *returnHTML = [NSMutableString string];

    int i = 0;
    NSError *error;
    NSDictionary *haplotypeJson = [NSJSONSerialization JSONObjectWithData:self.haplotypeData options:kNilOptions error:&error];
    NSUInteger haplotypeCount = haplotypeJson.count;
    [returnHTML appendString:@"<div class='datum-depth'>"];
    for (NSDictionary *haplotype in haplotypeJson) {
//        NSUInteger order = [[haplotype valueForKey:@"order"] unsignedIntegerValue];
        NSNumber *order = [hapReorder objectForKey:[haplotype valueForKey:@"order"]];
        NSString *depthString = [NSString stringWithFormat:@"%@", [haplotype valueForKey:@"depth"]];

        NSMutableString *appendString = [NSMutableString stringWithString:depthString];
        NSString *colorString = [self generateColorStringForOrder:order.unsignedIntegerValue];
        [appendString insertString:[NSString stringWithFormat:@"<font color='#%@'>", colorString] atIndex:0];
        [appendString appendString:@"</font>"];
//        NSDictionary *fontAttributes = [NSDictionary dictionaryWithObjectsAndKeys:[NSFont fontWithName:@"Courier" size:14], NSFontAttributeName, nil];

        [returnHTML appendString:appendString];
        if (i < haplotypeCount - 1) {
            [returnHTML appendString:@" / "];
        }
        ++i;
    }
    [returnHTML appendString:@"</div>"];

    return returnHTML;
}

- (NSMutableString *)renderNameHtml {
    NSString *formattedString = [self.name stringByReplacingOccurrencesOfString:@"_" withString:@" "];
    formattedString = [formattedString capitalizedString];

    NSMutableString *returnHTML = [NSMutableString stringWithFormat:@"<div class='datum-name'>%@</div>", formattedString];

    return returnHTML;
}

/**
* Gets an arbitrary allele array
*/
- (NSDictionary *)getHaplotypeOrder:(NSMutableArray *)alleleArray {
    NSError *error;
    NSMutableDictionary *haplotypeOrder = [NSMutableDictionary dictionaryWithCapacity:alleleArray.count];

    NSDictionary *haplotypeJson = [NSJSONSerialization JSONObjectWithData:self.haplotypeData options:kNilOptions error:&error];
    for (NSDictionary *haplotype in haplotypeJson) {
        NSString *haplotypeString = [NSString stringWithFormat:@"%@", [haplotype valueForKey:@"haplotype"]];
//        NSUInteger oldOrder = [[haplotype valueForKey:@"order"] unsignedIntegerValue];
        NSString* oldOrder = [haplotype valueForKey:@"order"] ;
        NSUInteger index =0 ;
        for (NSString *item in alleleArray) {
            if ([item rangeOfString:haplotypeString].location != NSNotFound) {
                [haplotypeOrder setObject:[NSNumber numberWithUnsignedInteger:index] forKey:oldOrder];
            }
            index++;
        }
//        NSString *colorString = [self generateColorStringForOrder:order];
    }

    return haplotypeOrder;
}
@end

