//
//  DatumMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "DatumMO.h"
#import "ConsensusStackEntryMO.h"
#import "DatumAlleleMO.h"
#import "DatumSnpMO.h"
#import "DepthMO.h"
#import "HaplotypeMO.h"
#import "LocusMO.h"
#import "ModelStackEntryMO.h"
#import "ReferenceStackEntryMO.h"
#import "SampleMO.h"
#import "StackEntryMO.h"


NSDictionary *getColorForOrder(NSUInteger order);

@implementation DatumMO

@dynamic name;
@dynamic sampleId;
@dynamic tagId;
@dynamic alleles;
@dynamic consensus;
@dynamic depths;
@dynamic haplotypes;
@dynamic locus;
@dynamic model;
@dynamic reference;
@dynamic sample;
@dynamic snps;
@dynamic stackEntries;


- (NSAttributedString *)renderName {
    NSString *formattedString = [self.name stringByReplacingOccurrencesOfString:@"_" withString:@" "];
    formattedString = [formattedString capitalizedString];

    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:formattedString];

    NSMutableParagraphStyle *mutParaStyle = [[NSMutableParagraphStyle alloc] init];
    [mutParaStyle setAlignment:NSCenterTextAlignment];
    [string addAttributes:[NSDictionary dictionaryWithObject:mutParaStyle
                                                      forKey:NSParagraphStyleAttributeName]
                    range:NSMakeRange(0, [[string string] length])];

    return string;
}

- (NSAttributedString *)renderHaplotypes {

    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] init];


//    for (NSUInteger i = 0; i < self.haplotypes.count; i++) {
    int i = 0;
//    NSArray *sortByOrder = [NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"order" ascending:YES]];

    NSArray *sortDescriptors = [NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"order" ascending:YES]];

    for (HaplotypeMO *haplotype  in [self.haplotypes sortedArrayUsingDescriptors:sortDescriptors]) {
//        NSString *haplotype = [self.haplotypes objectAtIndex:i];
        NSUInteger order = [self.locus lookupHaplotypeOrder:haplotype.haplotype];


        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:haplotype.haplotype];
        NSRange selectedRange = NSMakeRange(0, haplotype.haplotype.length);
        [appendString beginEditing];
        NSDictionary *attributes = getColorForOrder(order);
        [appendString setAttributes:attributes range:selectedRange];
        [appendString endEditing];

        [string appendAttributedString:appendString];
        if (i < self.haplotypes.count - 1) {
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

//    for(NSNumber *depth in self.depths){
//        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:[NSString stringWithFormat:@"%@",depth]];
//        [string appendAttributedString:appendString];
//    }
//    for (NSUInteger i = 0; i < self.depths.count; i++) {
    int i = 0;
    NSArray *sortDescriptors = [NSArray arrayWithObject:[NSSortDescriptor sortDescriptorWithKey:@"order" ascending:YES]];
    for (DepthMO *depth in  [self.depths sortedArrayUsingDescriptors:sortDescriptors]) {
//        NSString *depth = [(NSNumber *) [self.depths objectAtIndex:i] stringValue];

        NSUInteger order = [self.locus lookupDepthOrder:depth];


        NSString *numberString = [NSString stringWithFormat:@"%@", depth.depth];
        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:numberString];
        NSRange selectedRange = NSMakeRange(0, numberString.length);
        [appendString beginEditing];
        NSDictionary *attributes = getColorForOrder(order);
        [appendString setAttributes:attributes range:selectedRange];
        [appendString endEditing];

        [string appendAttributedString:appendString];
        if (i < self.depths.count - 1) {
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

- (NSArray *)getOrderedStackEntries {
    NSArray *sortedArray = [[self.stackEntries allObjects] sortedArrayUsingComparator:^NSComparisonResult(id a, id b) {
        StackEntryMO *first = (StackEntryMO*)a ;
        StackEntryMO *second = (StackEntryMO*)b ;
        NSString *firstRelationship = first.relationship ;
        NSString *secondRelationship = second.relationship ;

        int firstCount = 0 ;
        int secondCount = 0 ;
        //reference
        //consensus
        //model
        firstCount += [firstRelationship isEqualToString:@"reference"]?10000:0;
        secondCount += [secondRelationship isEqualToString:@"reference"]?10000:0;
        firstCount += [firstRelationship isEqualToString:@"consensus"]?1000:0;
        secondCount += [secondRelationship isEqualToString:@"consensus"]?1000:0;
        firstCount += [firstRelationship isEqualToString:@"model"]?100:0;
        secondCount += [secondRelationship isEqualToString:@"model"]?100:0;
        firstCount -= [first.entryId intValue];
        secondCount -= [second.entryId intValue];
        NSComparisonResult result = (firstCount>secondCount)?NSOrderedAscending:NSOrderedDescending;
        return result ;
    }];

    return sortedArray;

//    NSOrderedSet *orderedSet = [NSOrderedSet orderedSetWithArray:sortedArray];
//    return orderedSet;

}


@end

NSDictionary *getColorForOrder(NSUInteger order) {

    switch(order%3){
        case 0:
            return [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor greenColor], NSForegroundColorAttributeName,
                    nil];
        case 1:
            return [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor redColor], NSForegroundColorAttributeName,
                    nil];
        case 2:
            return [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor orangeColor], NSForegroundColorAttributeName,
                    nil];
        case 3:
            return [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor yellowColor], NSForegroundColorAttributeName,
                    nil];
        default:
            return [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor blueColor], NSForegroundColorAttributeName,
                    nil];
    }
}
