//
// Created by Nathan Dunn on 12/20/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//


#import "StackEntryRenderer.h"


@implementation StackEntryRenderer {

}

@synthesize sampleName;
@synthesize locusId;
@synthesize sequences;
@synthesize sequenceIds;
@synthesize blocks;
@synthesize relationships;
@synthesize entryIds;
@synthesize model;
@synthesize consensus;
@synthesize snpLocusData;
@synthesize snpDatumData;

@synthesize snpLocusLookup;
@synthesize snpDatumLookup;

- (id)init {
    self = [super init];
    if (self) {
        sequences = [[NSMutableArray alloc] init];
        sequenceIds = [[NSMutableArray alloc] init];
        blocks = [[NSMutableArray alloc] init];
        relationships = [[NSMutableArray alloc] init];
        entryIds = [[NSMutableArray alloc] init];
        snpLocusLookup = [[NSMutableDictionary alloc] init];
        snpDatumLookup = [[NSMutableDictionary alloc] init];
    }

    return self;
}


- (NSString *)renderHtml {

    NSError *error;
    NSDictionary *snpJson = [NSJSONSerialization JSONObjectWithData:snpLocusData options:kNilOptions error:&error];
    for (NSDictionary *snp in snpJson) {
        [snpLocusLookup setObject:snp forKey:[snp valueForKey:@"column"]];
    }

//    NSString* headerHTML= @"<head><link rel='stylesheet' type='text/css' href='test.css'/></head><body>" ;
//    NSString* cssPath = [[NSBundle mainBundle] pathForResource:@"test" ofType:@"css"];
    NSString *cssPath = [[NSBundle mainBundle] pathForResource:@"stacks" ofType:@"css"];
//    NSLog(@"css path %@",cssPath) ;
    NSString *cssString = [NSString stringWithContentsOfFile:cssPath encoding:NSUTF8StringEncoding error:NULL];
//    NSLog(@"css string %@",cssString) ;

//    NSString* returnHTML = [NSString stringWithFormat:@"<style type='text/css'>%@</style><div class='sample'>testy</div><table><tr><td col=4><div class='sample'>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</div></td></tr>",cssString,sampleName,locusId,[sequences count]];

//    NSMutableString* returnHTML = [NSMutableString stringWithFormat:@"<style type='text/css'>%@</style><div class='sample'>testy</div><table><tr><td col=4><div class='sample'>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</div></td></tr>",cssString,sampleName,locusId,[sequences count]];
    NSMutableString *returnHTML = [NSMutableString stringWithFormat:@"<style type='text/css'>%@</style><table class='radtag'>", cssString];
    [returnHTML appendString:[self renderHeader]];
    [returnHTML appendString:[self renderReference]];
    [returnHTML appendString:[self renderConsensus]];
    [returnHTML appendString:[self renderModel]];
    [returnHTML appendString:[self renderSequences]];
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


- (NSString *)renderReference {
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

- (NSString *)renderConsensus {
    NSMutableString *returnString = [NSMutableString string];
    NSMutableString *consensusString = [NSMutableString stringWithString:consensus];



    // TODO: sort by NSDictionary
    NSSortDescriptor *sortDescriptor = [NSSortDescriptor sortDescriptorWithKey:nil ascending:NO selector:@selector(localizedCompare:)];
    for (NSString *snpKey in [[snpLocusLookup allKeys] sortedArrayUsingDescriptors:[NSArray arrayWithObject:sortDescriptor]]) {
//        for (id snp in [snpJson allKeys]){
//        if(snpJson.count>1){
//            NSLog(@"snp %@ ",snpKey);
//        }
        NSUInteger column = [[NSNumber numberWithInteger:snpKey.integerValue] unsignedIntegerValue];
        [consensusString insertString:@"</span>" atIndex:(column + 1)];
        [consensusString insertString:@"<span class='rank_1'>" atIndex:(column)];
    }

    [returnString appendFormat:@"<tr><td class='num'></td><td class='con'>consensus</td><td class='id'></td><td class='tag'>%@</td></tr>", consensusString];
    return returnString;
}

- (NSString *)renderModel {
    NSMutableString *returnString = [NSMutableString string];
    [returnString appendFormat:@"<tr><td class='num'></td><td class='con'>model</td><td class='id'></td><td class='tag'>%@</td></tr>", model];
    return returnString;
}


- (NSMutableString *)renderSequences {
//    NSLog(@"sequences %ld entryIds %ld relatinships %ld sequenceIds %ld",sequences.count,entryIds.count,relationships.count,sequenceIds.count) ;

    NSMutableString *returnString = [NSMutableString string];

    NSString *blockStyle;
    for (NSUInteger  i = 0; i < sequences.count; i++) {

        // handle BLOCKS
        if([[blocks objectAtIndex:i] integerValue]%2 ==0 ){
            blockStyle =@"";
        }
        else{
            blockStyle =@" style='background-color: #dddddd;' ";
        }
        NSString* sequenceString = [sequences objectAtIndex:i] ;
        NSMutableString *formattedSequenceString = [NSMutableString stringWithString:sequenceString];

        // a dictionary of locations (integer) and formats to apply in reverse order
        NSMutableDictionary *formatDictionary = [NSMutableDictionary dictionary];

        // handle SNPS
        // for each locus snp, if a datum snp at the same location then


        // handle errors
        if([consensus isNotEqualTo:sequenceString]){
            for (int i = 0; i < sequenceString.length && i < consensus.length; i++) {
                if([ consensus characterAtIndex:i] != [sequenceString characterAtIndex:i]){
                    if([snpLocusLookup objectForKey:[NSNumber numberWithInt:i]]==nil){
                        [formatDictionary setObject:@"err" forKey:[NSNumber numberWithInt:i]];
                    }
                }
            }
        }
        
        // apply for the formats!!
        NSSortDescriptor *sortDescriptor = [NSSortDescriptor sortDescriptorWithKey:nil ascending:NO selector:@selector(compare:)];
        for (NSNumber *formatKey in [[formatDictionary allKeys] sortedArrayUsingDescriptors:[NSArray arrayWithObject:sortDescriptor]]) {
//            NSUInteger column = [[NSNumber numberWithInteger:formatKey.integerValue] unsignedIntegerValue];
            NSUInteger column = [formatKey unsignedIntegerValue];
            NSString* value = [formatDictionary objectForKey:formatKey];
            [formattedSequenceString insertString:@"</span>" atIndex:(column + 1)];
            [formattedSequenceString insertString:[NSString stringWithFormat:@"<span class='%@'>",value] atIndex:(column)];
        }



        [returnString appendFormat:@"<tr><td class='num'>%@</td><td class='%@'>%@</td><td class='id'>%@</td><td class='tag' %@>%@</td></tr>", [entryIds objectAtIndex:i], [relationships objectAtIndex:i], [relationships objectAtIndex:i], [sequenceIds objectAtIndex:i], blockStyle, formattedSequenceString];
    }
    return returnString;
}

- (BOOL)isEmpty {
    return sequences.count == 0;
}

- (void)clear {
    sampleName = @"";
    consensus = @"";
    model = @"";

    [sequences removeAllObjects];
    [sequenceIds removeAllObjects];
    [blocks removeAllObjects];
    [relationships removeAllObjects];
    [entryIds removeAllObjects];
    snpLocusData = nil ;
    [snpLocusLookup removeAllObjects];
    snpDatumData = nil ;
    [snpDatumLookup removeAllObjects];
}
@end