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
    for (HaplotypeMO *haplotype  in self.haplotypes) {
//        NSString *haplotype = [self.haplotypes objectAtIndex:i];
        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:haplotype.haplotype];
        NSRange selectedRange = NSMakeRange(0, haplotype.haplotype.length);
        [appendString beginEditing];
        NSDictionary *attributes;
        if (i % 2 == 0) {
            attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor greenColor], NSForegroundColorAttributeName,
                    nil];
        }
        else {
            attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor redColor], NSForegroundColorAttributeName,
                    nil];
        }
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
    for (DepthMO *depth in  self.depths) {
//        NSString *depth = [(NSNumber *) [self.depths objectAtIndex:i] stringValue];
        NSString *numberString = [NSString stringWithFormat:@"%@", depth.depth];
        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:numberString];
        NSRange selectedRange = NSMakeRange(0, numberString.length);
        [appendString beginEditing];
        NSDictionary *attributes;
        if (i % 2 == 0) {
            attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor greenColor], NSForegroundColorAttributeName,
                    nil];
        }
        else {
            attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSColor redColor], NSForegroundColorAttributeName,
                    nil];
        }
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
        //reference
        //consensus
        //model
        if([firstRelationship isEqualToString:@"reference"]){
            return NSOrderedAscending;
        }
        if([secondRelationship isEqualToString:@"reference"]){
            return NSOrderedDescending;
        }
        if([firstRelationship isEqualToString:@"consensus"]){
            return ([secondRelationship isEqualToString:@"reference"]?NSOrderedDescending:NSOrderedAscending);
        }
        if([firstRelationship isEqualToString:@"model"]){
            return ([secondRelationship isEqualToString:@"consensus"]||[secondRelationship isEqualToString:@"reference"]?NSOrderedDescending:NSOrderedAscending);
         }
        
//        if([secondRelationship isEqualToString:@"primary"] && firstRelationship
//           ){
//            return (![secondRelationship isEqualToString:@"model"]?NSOrderedAscending:NSOrderedDescending);
//        }
        // primary / secondary (just use entry id)
        NSComparisonResult result = [first.entryId compare:second.entryId];
        return result ;

    }];

    return sortedArray;

//    NSOrderedSet *orderedSet = [NSOrderedSet orderedSetWithArray:sortedArray];
//    return orderedSet;

}


@end
