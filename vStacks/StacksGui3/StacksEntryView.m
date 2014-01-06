//
// Created by Nathan Dunn on 12/20/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//


#import "StacksEntryView.h"


@implementation StacksEntryView {

}

@synthesize sampleName ;
@synthesize locusId;
@synthesize sequences;
@synthesize sequenceIds;
@synthesize blocks ;
@synthesize relationships;
@synthesize entryIds;

- (id)init {
    self = [super init];
    if (self) {
        sequences = [[NSMutableArray alloc] init];
        sequenceIds = [[NSMutableArray alloc] init];
        blocks = [[NSMutableArray alloc] init];
        relationships = [[NSMutableArray alloc] init];
        entryIds = [[NSMutableArray alloc] init];
    }

    return self;
}


- (NSString *)renderHtml {
//    NSString* headerHTML= @"<head><link rel='stylesheet' type='text/css' href='test.css'/></head><body>" ;
    NSString* cssPath = [[NSBundle mainBundle] pathForResource:@"test" ofType:@"css"];
//    NSLog(@"css path %@",cssPath) ;
    NSString *cssString = [NSString stringWithContentsOfFile:cssPath encoding:NSUTF8StringEncoding error:NULL];
//    NSLog(@"css string %@",cssString) ;

//    NSString* returnHTML = [NSString stringWithFormat:@"<style type='text/css'>%@</style><div class='sample'>testy</div><table><tr><td col=4><div class='sample'>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</div></td></tr>",cssString,sampleName,locusId,[sequences count]];
    
//    NSMutableString* returnHTML = [NSMutableString stringWithFormat:@"<style type='text/css'>%@</style><div class='sample'>testy</div><table><tr><td col=4><div class='sample'>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</div></td></tr>",cssString,sampleName,locusId,[sequences count]];
    NSMutableString* returnHTML = [NSMutableString stringWithFormat:@"<style type='text/css'>%@</style><table>",cssString];
    [returnHTML appendString:[self renderReference]];
    [returnHTML appendString:[self renderConsensus]];
    [returnHTML appendString:[self renderModel]];
    [returnHTML appendString:[self renderSequences]];
    [returnHTML appendString:@"</table></body>"];
    return returnHTML ;
}


- (NSString *)renderReference {
    return @"";
}

- (NSString *)renderConsensus {
    return @"";
}

- (NSString *)renderModel {
    return @"";
}


- (NSMutableString *)renderSequences {
    NSMutableString* returnString = [[NSMutableString alloc] init];
//    NSString *sequence  ;
//    NSString *sequenceId  ;
    for(int i = 0 ; i < sequences.count ; i++){
        [returnString appendFormat:@"<tr><td>%@</td><td>%@</td><td>%@</td><td>%@</td></tr>",[entryIds objectAtIndex:i],[relationships objectAtIndex:i],[sequenceIds objectAtIndex:i],[sequences objectAtIndex:i]];
    }
    return returnString;
}

- (BOOL)isEmpty {
    return sequences.count==0;
}

- (void)clear {
    [sequences removeAllObjects];
    [sequenceIds removeAllObjects];
}
@end