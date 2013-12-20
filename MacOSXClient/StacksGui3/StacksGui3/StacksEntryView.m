//
// Created by Nathan Dunn on 12/20/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//


#import "StacksEntryView.h"


@implementation StacksEntryView {

}

@synthesize sampleName ;
@synthesize locusId;
@synthesize sequenceIds;
@synthesize sequences;

- (id)init {
    self = [super init];
    if (self) {
        sequences = [[NSMutableArray alloc] init];
        sequenceIds = [[NSMutableArray alloc] init];
    }

    return self;
}


- (NSString *)renderHtml {
     NSString* returnHTML = [NSString stringWithFormat:@"<table><tr><td col=4>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</td></tr>",sampleName,locusId,[sequences count]];

//    returnHTML = [NSString stringWithFormat:@"%@%@",returnHTML,[self renderSequences]]; ;


    NSString* sequenceString = [self renderSequences];

    return [NSString stringWithFormat:@"%@%@</table>",returnHTML,sequenceString];
}

- (NSString *)renderSequences {
    NSMutableString* returnString = [NSMutableString string];

//    NSString *sequence ;
//    NSString *sequenceId ;
    for(int i = 0 ; i < sequences.count  ; i++){
//        sequence = [sequences objectAtIndex:i];
//        sequenceId = [sequenceIds  objectAtIndex:i];
//        returnString = [NSString stringWithFormat:@"%@<tr><td>%@</td><td col=3>%@</td></tr>",returnString,[sequenceIds objectAtIndex:i]
//                ,[sequences objectAtIndex:i]];
        [returnString appendFormat:@"<tr><td>%@</td><td col=3>%@</td></tr>",[sequenceIds objectAtIndex:i]
                ,[sequences objectAtIndex:i]];
    }
    return returnString;
}
@end