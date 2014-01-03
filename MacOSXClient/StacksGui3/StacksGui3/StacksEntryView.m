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

- (id)init {
    self = [super init];
    if (self) {
        sequences = [[NSMutableArray alloc] init];
        sequenceIds = [[NSMutableArray alloc] init];
    }

    return self;
}



- (NSString *)renderHtml {
//    NSString* headerHTML= @"<head><link rel='stylesheet' type='text/css' href='test.css'/></head><body>" ;
//    NSString* cssPath = [[NSBundle mainBundle] pathForResource:@"test" ofType:@"css"];
//    NSLog(@"css path %@",cssPath) ;
//    NSString *cssString = [NSString stringWithContentsOfFile:cssPath encoding:NSUTF8StringEncoding error:NULL];
    
//    NSLog(@"css string %@",cssString) ;

//    NSString* returnHTML = [NSString stringWithFormat:@"<style type='text/css'>%@</style><div class='sample'>testy</div><table><tr><td col=4><div class='sample'>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</div></td></tr>",cssString,sampleName,locusId,[sequences count]];

    NSMutableString* returnHTML = [NSMutableString stringWithString:@"<style type='text/css'>"];
    [returnHTML appendString:self.cssString];
    [returnHTML appendString:@"</style><div class='sample'>testy</div><table><tr><td col=4><div class='sample'>RENDERED Some stack data for sample "];
    [returnHTML appendString:sampleName];
    [returnHTML appendString:@"</div></td></tr>"];

//    [returnHTML appendString:@"</style><div class='sample'>testy</div><table><tr><td col=4><div class='sample'>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</div></td></tr>"];
//    NSMutableString* returnHTML = [NSMutableString stringWithFormat:@"<style type='text/css'>%@</style><div class='sample'>testy</div><table><tr><td col=4><div class='sample'>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</div></td></tr>",cssString,sampleName,locusId,[sequences count]];

//    returnHTML = [NSString stringWithFormat:@"%@%@",returnHTML,[self renderSequences]]; ;
    [returnHTML appendString:[self renderSequences]];
    [returnHTML appendString:@"</table></body>"];
    return returnHTML ;
}

- (NSMutableString *)renderSequences {
    NSMutableString* returnString = [[NSMutableString alloc] init];
//    NSString *sequence  ;
//    NSString *sequenceId  ;
    for(int i = 0 ; i < sequences.count ; i++){
        [returnString appendFormat:@"<tr><td>%@</td><td col=3>%@</td></tr>",[sequences objectAtIndex:i],[sequenceIds objectAtIndex:i]];
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