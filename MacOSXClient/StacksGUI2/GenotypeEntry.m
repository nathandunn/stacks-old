//
// Created by ndunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "GenotypeEntry.h"


@implementation GenotypeEntry {

}

//@synthesize type;
//@synthesize superScript;
//@synthesize subScript;


- (NSString *)render {
    NSString *renderString = @"";
    for(int i = 0 ; i < _haplotypes.count ; i++){
        renderString = [renderString stringByAppendingString:[_haplotypes objectAtIndex:i]];
    }
    for(int i = 0 ; i < _depths.count ; i++){
        renderString = [renderString stringByAppendingString:[_depths objectAtIndex:i]];
    }
//    NSLog(@"rendering string %@",renderString);


    return renderString;
}
@end