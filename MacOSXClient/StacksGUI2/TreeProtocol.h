//
// Created by ndunn on 4/1/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@protocol TreeProtocol <NSObject>


- (NSUInteger)childCount;
- (id)childAtIndex:(NSUInteger) index;
- (BOOL)isLeaf;

@end
