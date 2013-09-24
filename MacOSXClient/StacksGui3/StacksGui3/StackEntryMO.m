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


@implementation StackEntryMO

@dynamic block;
@dynamic entryId;
@dynamic relationship;
@dynamic sequence;
@dynamic sequenceId;
@dynamic datum;


- (NSAttributedString *)renderEntryId{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:[NSString stringWithFormat:@"%@",self.entryId] ];
    return string ;
}

- (NSAttributedString *)renderRelationship{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:[NSString stringWithFormat:@"%@",self.relationship] ];
    [string beginEditing];
    NSDictionary *attributes ;

    if([self.relationship isEqualToString:@"primary"]){
    attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor greenColor], NSForegroundColorAttributeName,
            nil];
    }
    // secondary
    else{
        attributes = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor redColor], NSForegroundColorAttributeName,
                nil];
    }
    [string setAttributes:attributes range:NSMakeRange(0, self.relationship.length)];
    [string endEditing];
    return string ;
}

- (NSAttributedString *)renderSequence {
//    NSMutableAttributedString *attributedString = [[NSMutableAttributedString alloc] initWithString:self.sequence];

    NSMutableAttributedString *sequenceAttributedString = [[NSMutableAttributedString alloc] initWithString:self.sequence];

    [sequenceAttributedString beginEditing];
    NSDictionary *blockAttribute;
    NSDictionary *snpAttribute;
    NSDictionary *defectAttribute;


    if(![self.block isEqualToString:@"1"]){
        snpAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blackColor], NSForegroundColorAttributeName,
                [NSColor controlShadowColor], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
                nil];
        blockAttribute= [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blackColor], NSForegroundColorAttributeName,
                [NSColor whiteColor], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
                nil];
        defectAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blackColor], NSForegroundColorAttributeName,
                [NSColor redColor], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
                nil];
    }
    else{
        snpAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor redColor], NSForegroundColorAttributeName,
                [NSColor lightGrayColor], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
                nil];
        blockAttribute= [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor blackColor], NSForegroundColorAttributeName,
                [NSColor grayColor], NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
                nil];
        defectAttribute = [NSDictionary dictionaryWithObjectsAndKeys:
                [NSColor whiteColor], NSForegroundColorAttributeName,
                [NSColor redColor] , NSBackgroundColorAttributeName,
                [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
                nil];
    }

    [sequenceAttributedString setAttributes:blockAttribute range:NSMakeRange(0, self.sequence.length)];

    NSMutableArray* snps = [[NSMutableArray alloc] init];
    for (LocusSnpMO *snp in self.datum.locus.snps) {
        NSRange selectedRange = NSMakeRange([snp.column unsignedIntegerValue], 1);
        [sequenceAttributedString setAttributes:snpAttribute range:selectedRange];
        [snps addObject:snp.column];
    }

    NSString* consensusSequence = self.datum.consensus.sequence;
    NSString* mySequence = [sequenceAttributedString string];
    if(![consensusSequence isEqualToString:mySequence]){
        for(int i = 0 ; i < mySequence.length && i < consensusSequence.length ; i++){
            if([consensusSequence characterAtIndex:i]!=[mySequence characterAtIndex:i]){
                // if is actually a defined SNP, then we ignore
                if(![snps containsObject:[NSNumber numberWithInt:i]]){
                    NSRange selectedRange = NSMakeRange(i, 1);
                    [sequenceAttributedString setAttributes:defectAttribute range:selectedRange];
                }
            }
        }
    }



    [sequenceAttributedString endEditing];
    return sequenceAttributedString;


//    return attributedString ;
}
@end
