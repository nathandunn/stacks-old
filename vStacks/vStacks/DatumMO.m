//
//  DatumMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
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
#import "SampleMO.h"
#import "PopulationMO.h"


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

@synthesize colorGenerator;


- (NSAttributedString *)renderName {

    PopulationMO *populationMO = self.sample.population;
    NSString* populationName = @"";
    if(populationMO!=nil) {
        populationName = [populationName stringByAppendingFormat:@": %@",populationName];
    }
    NSLog(@"population name: %@",populationName) ;


    NSString *formattedString = [self.name stringByReplacingOccurrencesOfString:@"_" withString:@" "];
    formattedString = [formattedString capitalizedString];
    formattedString = [formattedString stringByAppendingString:populationName];

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

- (NSAttributedString *)renderPopulation{
    PopulationMO *populationMO = self.sample.population;
    if(populationMO==nil) {
        return [[NSAttributedString alloc] init];
    }
    NSString *formattedString = [populationMO.name stringByReplacingOccurrencesOfString:@"_" withString:@" "];
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


//    for (NSUInteger i = 0; i < self.haplotypes.count; i++) {
    int i = 0;
//    NSArray *sortByOrder = [NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"order" ascending:YES]];

//    NSArray *sortDescriptors = [NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"order" ascending:YES]];

    // TODO: convert
    NSError *error;
    NSDictionary *haplotypeJson = [NSJSONSerialization JSONObjectWithData:self.haplotypeData options:kNilOptions error:&error];
    NSUInteger haplotypeCount = haplotypeJson.count;
    for (NSDictionary *haplotype in haplotypeJson) {
//    for (HaplotypeMO *haplotype  in [self.haplotypes sortedArrayUsingDescriptors:sortDescriptors]) {
////        NSString *haplotype = [self.haplotypes objectAtIndex:i];

//        NSUInteger order = [self.locus lookupHaplotypeOrder:haplotype.haplotype];
        NSUInteger order = [[haplotype valueForKey:@"order"] unsignedIntegerValue];
        NSString *haplotypeString = [haplotype valueForKey:@"haplotype"];
//        NSUInteger depth = [[haplotype valueForKey:@"depth"] unsignedIntegerValue];


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
        if (i < haplotypeCount - 1 ) {
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
    NSError *error ;
    NSDictionary *haplotypeJson = [NSJSONSerialization JSONObjectWithData:self.haplotypeData options:kNilOptions error:&error];
    NSUInteger haplotypeCount = haplotypeJson.count;
    for (NSDictionary *haplotype in haplotypeJson) {
//    for (HaplotypeMO *haplotype  in [self.haplotypes sortedArrayUsingDescriptors:sortDescriptors]) {
////        NSString *haplotype = [self.haplotypes objectAtIndex:i];

//        NSUInteger order = [self.locus lookupHaplotypeOrder:haplotype.haplotype];
        NSUInteger order = [[haplotype valueForKey:@"order"] unsignedIntegerValue];
//        NSString *haplotypeString = [haplotype valueForKey:@"haplotype"];
//        NSUInteger depth = [[haplotype valueForKey:@"depth"] unsignedIntegerValue];
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
//        if (i < self.haplotypes.count - 1) {
        if (i < haplotypeCount - 1 ) {
            [string appendAttributedString:[[NSAttributedString alloc] initWithString:@" / "]];
        }
        ++i;
    }



//    NSArray *sortDescriptors = [NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"order" ascending:YES]];
    // TODO: convert
//    for (DepthMO *depth in  [self.depths sortedArrayUsingDescriptors:sortDescriptors]) {
////        NSString *depth = [(NSNumber *) [self.depths objectAtIndex:i] stringValue];
//
//        NSUInteger order = [self.locus lookupDepthOrder:depth];
//
//
//        NSString *numberString = [NSString stringWithFormat:@"%@", depth.depth];
//        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:numberString];
//        NSRange selectedRange = NSMakeRange(0, numberString.length);
//        [appendString beginEditing];
//        NSDictionary *attributes = [self generateColorForOrder:order];
//        [appendString setAttributes:attributes range:selectedRange];
//        NSDictionary *fontAttributes = [NSDictionary dictionaryWithObjectsAndKeys:[NSFont fontWithName:@"Courier" size:14], NSFontAttributeName, nil];
//        [appendString addAttributes:fontAttributes range:selectedRange];
//        [appendString endEditing];
//
//        [string appendAttributedString:appendString];
//        if (i < self.depths.count - 1) {
//            [string appendAttributedString:[[NSAttributedString alloc] initWithString:@" / "]];
//        }
//        ++i;
//    }


    NSMutableParagraphStyle *mutParaStyle = [[NSMutableParagraphStyle alloc] init];
    [mutParaStyle setAlignment:NSCenterTextAlignment];
    [string addAttributes:[NSDictionary dictionaryWithObject:mutParaStyle
                                                      forKey:NSParagraphStyleAttributeName]
                    range:NSMakeRange(0, [[string string] length])];

    return string;
}

//- (NSArray *)getOrderedStackEntries {
//    NSArray *sortedArray = [[self.stackEntries allObjects] sortedArrayUsingComparator:^NSComparisonResult(id a, id b) {
//        StackEntryMO *first = (StackEntryMO *) a;
//        StackEntryMO *second = (StackEntryMO *) b;
//        NSString *firstRelationship = first.relationship;
//        NSString *secondRelationship = second.relationship;
//
//        int firstCount = 0;
//        int secondCount = 0;
//        //reference
//        //consensus
//        //model
//        firstCount += [firstRelationship isEqualToString:@"reference"] ? 10000 : 0;
//        secondCount += [secondRelationship isEqualToString:@"reference"] ? 10000 : 0;
//        firstCount += [firstRelationship isEqualToString:@"consensus"] ? 1000 : 0;
//        secondCount += [secondRelationship isEqualToString:@"consensus"] ? 1000 : 0;
//        firstCount += [firstRelationship isEqualToString:@"model"] ? 100 : 0;
//        secondCount += [secondRelationship isEqualToString:@"model"] ? 100 : 0;
//        firstCount -= [first.entryId intValue];
//        secondCount -= [second.entryId intValue];
//        NSComparisonResult result = (firstCount > secondCount) ? NSOrderedAscending : NSOrderedDescending;
//        return result;
//    }];
//
//    return sortedArray;
//
////    NSOrderedSet *orderedSet = [NSOrderedSet orderedSetWithArray:sortedArray];
////    return orderedSet;
//
//}

- (NSDictionary *)generateColorForOrder:(NSUInteger)order {
    if (colorGenerator == nil) {
        colorGenerator = [[ColorGenerator alloc] init];
    }

    NSColor *color = [colorGenerator generateColorForOrder:order];
    return [NSDictionary dictionaryWithObjectsAndKeys:
            color, NSForegroundColorAttributeName,
            nil];

}

@end

