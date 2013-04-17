//
// Created by ndunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "GenotypeEntry.h"


@implementation GenotypeEntry {

}

//@synthesize tagId;
//@synthesize superScript;
//@synthesize subScript;

@synthesize name = _name;

- (id)init {
    self = [super init];
    if (self) {
        _name = @"Test";
    }
    return self;
}


- (id)initWithCapcity:(int)numHaplotypes {
    self = [super init];
    if (self) {
        _name = @"Test";
        _haplotypes = [[NSMutableArray alloc] initWithCapacity:numHaplotypes];
        _depths = [[NSMutableArray alloc] initWithCapacity:numHaplotypes];
    }
    return self;
}


- (NSString *)render {
    NSString *renderString = @"";
    for (int i = 0; i < _haplotypes.count; i++) {
        renderString = [renderString stringByAppendingString:[_haplotypes objectAtIndex:i]];
    }
    for (int i = 0; i < _depths.count; i++) {
        renderString = [renderString stringByAppendingString:[_depths objectAtIndex:i]];
    }
//    NSLog(@"rendering string %@",renderString);


    return renderString;
}


- (NSAttributedString *)renderName {
    NSString *formattedString = [_name stringByReplacingOccurrencesOfString:@"_" withString:@" "];
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

    for (NSUInteger i = 0; i < self.haplotypes.count; i++) {
        NSString *haplotype = [self.haplotypes objectAtIndex:i];
        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:haplotype];
        NSRange selectedRange = NSMakeRange(0, haplotype.length);
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
    for (NSUInteger i = 0; i < self.depths.count; i++) {
        NSString *depth = [(NSNumber *) [self.depths objectAtIndex:i] stringValue];
        NSMutableAttributedString *appendString = [[NSMutableAttributedString alloc] initWithString:depth];
        NSRange selectedRange = NSMakeRange(0, depth.length);
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
    }

    NSMutableParagraphStyle *mutParaStyle = [[NSMutableParagraphStyle alloc] init];
    [mutParaStyle setAlignment:NSCenterTextAlignment];
    [string addAttributes:[NSDictionary dictionaryWithObject:mutParaStyle
                                                      forKey:NSParagraphStyleAttributeName]
                    range:NSMakeRange(0, [[string string] length])];

    return string;
}

@end

