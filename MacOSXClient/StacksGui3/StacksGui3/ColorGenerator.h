//
//  ColorGenerator.h
//  StacksGui3
//
//  Created by Nathan Dunn on 10/22/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface ColorGenerator : NSObject

- (NSColor *)generateColorForOrder:(NSUInteger)order;
@end
