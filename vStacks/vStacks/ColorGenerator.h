//
//  ColorGenerator.h
//  StacksGui3
//
//  Created by Nathan Dunn on 10/22/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface ColorGenerator : NSObject

- (NSColor *)generateColorForOrder:(NSUInteger)order;
- (NSColor *)colorWithHexString:(NSString *)hex;
- (NSString *)generateColorStringForOrder:(NSUInteger)order ;

@end
