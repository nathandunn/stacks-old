//
// Created by Nathan Dunn on 1/14/14.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//


#import "StacksEntryDatumRenderer.h"


@implementation StacksEntryDatumRenderer {

}

@synthesize snpLocusLookup;
@synthesize snpDatumLookup;

@synthesize numberFormatter;

- (id)init {
    self = [super init];
    if (self) {
        snpLocusLookup = [NSMutableDictionary dictionary];
        snpDatumLookup = [NSMutableDictionary dictionary];


        numberFormatter = [[NSNumberFormatter alloc] init];
        numberFormatter.numberStyle = NSNumberFormatterNoStyle;
    }

    return self;
}


- (NSString *)renderHtmlForData:(NSData *)data datumSnps:(NSData *)snpDatumData locusSnps:(NSData *)snpLocusData {

    NSError *error;
    if (snpLocusData != nil ) {
        NSDictionary *snpJson = [NSJSONSerialization JSONObjectWithData:snpLocusData options:kNilOptions error:&error];
        for (NSDictionary *snp in snpJson) {
            // not sure why this is a string, clearly an NSNumber when put in . . . weird
            [snpLocusLookup setObject:snp forKey:[numberFormatter numberFromString:[snp valueForKey:@"column"]]];
//            [snpLocusLookup setObject:snp forKey:[snp valueForKey:@"column"]];
        }

    }

    if (snpDatumData != nil) {
        NSDictionary *snpJson = [NSJSONSerialization JSONObjectWithData:snpDatumData options:kNilOptions error:&error];
        for (NSDictionary *snp in snpJson) {
            [snpDatumLookup setObject:snp forKey:[snp valueForKey:@"column"]];
        }
    }

    NSString *cssPath = [[NSBundle mainBundle] pathForResource:@"stacks" ofType:@"css"];
    NSString *cssString = [NSString stringWithContentsOfFile:cssPath encoding:NSUTF8StringEncoding error:NULL];

    NSMutableString *returnHTML = [NSMutableString stringWithFormat:@"<style type='text/css'>%@</style><table class='radtag'>", cssString];


    NSDictionary *jsonData = [NSJSONSerialization JSONObjectWithData:data options:kNilOptions error:&error];
    NSString *consensus = [[jsonData objectForKey:@"consensus"] objectForKey:@"sequence"];
    NSString *model = [[jsonData objectForKey:@"model"] objectForKey:@"sequence"];


    [returnHTML appendString:[self renderReference:consensus]];
    [returnHTML appendString:[self renderConsensus:consensus]];
    [returnHTML appendString:[self renderModel:model]];
    [returnHTML appendString:[self renderSequences:jsonData]];
    [returnHTML appendString:@"</table></body>"];
    return returnHTML;
}

- (NSString *)renderHeader {
    NSMutableString *returnString = [NSMutableString string];
    [returnString appendString:@"<tr>"];
    [returnString appendString:@"<th style='width: 5%;'>&nbsp;</th>"];
    [returnString appendString:@"<th style='width: 15%;'>Relationship</th>"];
    [returnString appendString:@"<th style='width: 20%;'>Seq ID</th>"];
    [returnString appendString:@"<th style='width: 60%;'>Sequence</th>"];
    [returnString appendString:@"</tr>"];
    return returnString;
}

- (NSString *)renderReference:(NSString *)consensus {
    NSMutableString *returnString = [NSMutableString string];
    NSMutableString *referenceString = [NSMutableString string];
    NSUInteger sequenceSize = consensus.length;

    BOOL chunk1 = 0;
    for (NSUInteger i = 0; i < sequenceSize; i++) {
        if (i % 10 == 0) {
            if (i > 0) {
                [referenceString appendString:@"</span>"];
            }
            chunk1 = !chunk1; // flip it
            [referenceString appendFormat:@"<span class='%@'>", (chunk1 ? @"light_scale" : @"dark_scale")];
        }
        [referenceString appendFormat:@"%ld", i % 10];
    }
    [referenceString appendString:@"</span>"];
    [returnString appendFormat:@"<tr><td class='num'></td><td class='con'></td><td class='id'></td><td class='tag'>%@</td></tr>", referenceString];
    return returnString;
}

- (NSString *)renderConsensus:(NSString *)consensus {
    NSMutableString *returnString = [NSMutableString string];
    NSMutableString *consensusString = [NSMutableString stringWithString:consensus];


    NSSortDescriptor *sortDescriptor = [NSSortDescriptor sortDescriptorWithKey:nil ascending:NO selector:@selector(compare:)];
    for (NSNumber *snpKey in [[snpLocusLookup allKeys] sortedArrayUsingDescriptors:[NSArray arrayWithObject:sortDescriptor]]) {
        NSUInteger column = snpKey.unsignedIntegerValue;
        if (column > consensusString.length) {
            NSLog(@"Consensus snp column %ld is greater %ld ", column, consensusString.length);
        }
        else {
            if (column >= consensusString.length - 1) {
                [consensusString appendString:@"</span>"];
            }
            else {
                [consensusString insertString:@"</span>" atIndex:(column + 1)];
            }
            [consensusString insertString:@"<span class='rank_1'>" atIndex:(column)];
        }
    }

    [returnString appendFormat:@"<tr><td class='num'></td><td class='con'>consensus</td><td class='id'></td><td class='tag'>%@</td></tr>", consensusString];
    return returnString;
}

- (NSString *)renderModel:(NSString *)model {
    NSMutableString *returnString = [NSMutableString string];
    [returnString appendFormat:@"<tr><td class='num'></td><td class='con'>model</td><td class='id'></td><td class='tag'>%@</td></tr>", model];
    return returnString;
}

- (NSMutableString *)renderSequences:(NSDictionary *)sequenceDictionary {
//    NSLog(@"sequences %ld entryIds %ld relatinships %ld sequenceIds %ld",sequences.count,entryIds.count,relationships.count,sequenceIds.count) ;

    NSString *consensus = [[sequenceDictionary objectForKey:@"consensus"] objectForKey:@"sequence"];

    NSUInteger numSequences = sequenceDictionary.count - 2;

    NSMutableString *returnString = [NSMutableString string];

    NSString *blockStyle;
    for (NSUInteger i = 0; i < numSequences; i++) {

        NSString *key = [NSNumber numberWithInt:i+1].stringValue;
        NSDictionary *sequenceObject = [sequenceDictionary objectForKey:key];
        NSString *block = [sequenceObject objectForKey:@"block"];

        // handle BLOCKS
        if ([block integerValue] % 2 == 0) {
            blockStyle = @"";
        }
        else {
            blockStyle = @" style='background-color: #dddddd;' ";
        }

        NSString *sequenceString = [sequenceObject objectForKey:@"sequence"];

        NSMutableString *formattedSequenceString = [NSMutableString stringWithString:sequenceString];

        // a dictionary of locations (integer) and formats to apply in reverse order
        NSMutableDictionary *formatDictionary = [NSMutableDictionary dictionary];
        NSSortDescriptor *sortDescriptor = [NSSortDescriptor sortDescriptorWithKey:nil ascending:NO selector:@selector(compare:)];

        // handle SNPS
        // for each locus snp, if a datum snp at the same location then

        for (NSNumber *snpColumn in [[snpLocusLookup allKeys] sortedArrayUsingDescriptors:[NSArray arrayWithObject:sortDescriptor]]) {
            NSUInteger column = snpColumn.unsignedIntegerValue;

            if ([snpDatumLookup objectForKey:snpColumn] != nil) {
                if ([consensus characterAtIndex:column] == [sequenceString characterAtIndex:column]) {
                    [formatDictionary setObject:@"rank_1" forKey:snpColumn];
                }
                else {
                    [formatDictionary setObject:@"rank_2" forKey:snpColumn];
                }
            }
            else {
                [formatDictionary setObject:@"cat_snp" forKey:snpColumn];
            }

        }


        // handle errors
        if ([consensus isNotEqualTo:sequenceString]) {
            for (int i = 0; i < sequenceString.length && i < consensus.length; i++) {
                if ([consensus characterAtIndex:i] != [sequenceString characterAtIndex:i]) {
                    if ([snpLocusLookup objectForKey:[NSNumber numberWithInt:i]] == nil && [formatDictionary objectForKey:[NSNumber numberWithInt:i]] == nil) {
                        [formatDictionary setObject:@"err" forKey:[NSNumber numberWithInt:i]];
                    }
                }
            }
        }

        // apply for the formats!!
        for (NSNumber *formatKey in [[formatDictionary allKeys] sortedArrayUsingDescriptors:[NSArray arrayWithObject:sortDescriptor]]) {
            NSUInteger column = [formatKey unsignedIntegerValue];
            NSString *value = [formatDictionary objectForKey:formatKey];

            if (column > formattedSequenceString.length) {
                NSLog(@"Bad format on sequence string %ld > %ld", column, formattedSequenceString.length);
            }
            else {
                if (column >= formattedSequenceString.length - 1) {
                    [formattedSequenceString appendString:@"</span>"];
                }
                else {
                    [formattedSequenceString insertString:@"</span>" atIndex:(column + 1)];
                }
                [formattedSequenceString insertString:[NSString stringWithFormat:@"<span class='%@'>", value] atIndex:(column)];
            }
        }


        [returnString appendFormat:@"<tr><td class='num'>%@</td><td class='%@'>%@</td><td class='id'>%@</td><td class='tag' %@>%@</td></tr>"
                , key
                , [sequenceObject objectForKey:@"relationship"]
                , [sequenceObject objectForKey:@"relationship"]
                , [sequenceObject objectForKey:@"sequenceId"]
                , blockStyle, formattedSequenceString];
    }
    return returnString;
}

@end