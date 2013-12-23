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
    }

    return self;
}


- (NSString *)renderHtml {
    NSString* returnHTML = [NSString stringWithFormat:@"<script type='text/javascript'></script> <table><tr><td col=4><div style='color:red;display:inline;' class='sample'>RENDERED Some stack data for sample '%@' and locus '%ld' and # stacks: %ld</div></td></tr>",sampleName,locusId,[sequences count]];

    returnHTML = [NSString stringWithFormat:@"%@%@",returnHTML,[self renderSequences]]; ;



    returnHTML = [NSString stringWithFormat:@"%@</table>",returnHTML];
    return returnHTML ;
}

- (NSString *)renderSequences {
    NSMutableString* returnString = [[NSMutableString alloc] init];
    NSString *sequence  ;
    NSString *sequenceId  ;
    for(int i = 0 ; i < sequences.count ; i++){
        sequence = [sequences objectAtIndex:i];
        sequenceId = [sequenceIds objectAtIndex:i];
        [returnString appendFormat:@"<tr><td col=4 style='font-family:monospace';>%@ - %@</td></tr>",sequenceId,sequence];
    }
//    for(NSString* sequence in sequences){
//    }
    return returnString;
}
@end