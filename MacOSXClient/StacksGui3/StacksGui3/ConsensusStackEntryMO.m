//
//  ConsensusStackEntryMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "ConsensusStackEntryMO.h"
#import "DatumSnpMO.h"
#import "DatumMO.h"


@implementation ConsensusStackEntryMO

- (NSAttributedString *)renderEntryId{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:@""];
    return string ;
}

- (NSAttributedString *)renderRelationship{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.relationship];
    return string ;
}

- (NSAttributedString *)renderSequence {
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.sequence];
    [string beginEditing];
    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor redColor], NSForegroundColorAttributeName,
            [NSColor grayColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier Bold" size:14.0], NSFontAttributeName,
            nil];
    for (DatumSnpMO *snp in self.datum.snps) {
        NSRange selectedRange = NSMakeRange([snp.column unsignedIntegerValue], 1);
        [string setAttributes:attributes range:selectedRange];
    }
    [string endEditing];
    return string;
}

@end
