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
    if(_superScript!=nil && _subScript!=nil){
        return [NSString stringWithFormat:@"%@ / %@",_superScript,_subScript];
    }
    else
    if(_superScript!=nil){
        return [NSString stringWithFormat:@"%@",_superScript];
    }
    else
    if(_subScript!=nil){
        return [NSString stringWithFormat:@"%@",_subScript];
    }

    return @"";
}
@end