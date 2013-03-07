//
// Created by ndunn on 2/27/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StacksView.h"
#import "StackEntry.h"


@implementation StacksView {


}

- (NSInteger)rowsNeeded {
    NSLog(@"rows needed calc") ;
 return 3+ [_stackEntries count];
}
@end