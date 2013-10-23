//
//  ColorGenerator.m
//  StacksGui3
//
//  Created by Nathan Dunn on 10/22/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
#define SRGB (CGFloat [4]){184.0/256.0, 184.0/256.0, 184.0/256.0, 1.0}


#import "ColorGenerator.h"

@implementation ColorGenerator

/*
#008000
#c00000
#ffc600
#29356c
#860000
#dc6200
#4b398e
#008f56
#bf1e25
#4cb8ff
 */
- (NSColor *)generateColorForOrder:(NSUInteger)order {
    switch (order % 10) {
        case 0:
            return [self colorWithHexString:@"008000"];
        case 1:
            return [self colorWithHexString:@"c00000"];
        case 2:
            return [self colorWithHexString:@"ffc600"];
        case 3:
            return [self colorWithHexString:@"29356c"];
        case 4:
            return [self colorWithHexString:@"860000"];
        case 5:
            return [self colorWithHexString:@"dc6200"];
        case 6:
            return [self colorWithHexString:@"4b398e"];
        case 7:
            return [self colorWithHexString:@"008f56"];
        case 8:
            return [self colorWithHexString:@"bf1e25"];
        case 9:
            return [self colorWithHexString:@"4cb8ff"];
        default:
            NSLog(@"should never be hear");
            return nil;
    }
}

- (NSColor *)colorWithHexStringA:(NSString *)hexString {
    unsigned rgbValue = 0;
    NSScanner *scanner = [NSScanner scannerWithString:hexString];
    [scanner setScanLocation:1]; // bypass '#' character
    [scanner scanHexInt:&rgbValue];
    return [NSColor colorWithDeviceRed:((rgbValue & 0xFF0000) >> 16)/255.0 green:((rgbValue & 0xFF00) >> 8)/255.0 blue:(rgbValue & 0xFF)/255.0 alpha:1.0];
}

- (NSColor *)colorWithHexString:(NSString *)hex {
//    [NSColorSpace genericRGBColorSpace];
    NSLog(@"hex %@", hex);
    NSString *cString = [[hex stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]] uppercaseString];
    NSLog(@"hex->cString %@->%@", hex, cString);

    // String should be 6 or 8 characters
    if ([cString length] < 6) return [NSColor grayColor];

    // strip 0X if it appears
    if ([cString hasPrefix:@"0X"]) cString = [cString substringFromIndex:2];

    if ([cString length] != 6) return [NSColor grayColor];

    // Separate into r, g, b substrings
    NSRange range;
    range.location = 0;
    range.length = 2;
    NSString *rString = [cString substringWithRange:range];

    range.location = 2;
    NSString *gString = [cString substringWithRange:range];

    range.location = 4;
    NSString *bString = [cString substringWithRange:range];

    // Scan values
    unsigned int r, g, b;
    [[NSScanner scannerWithString:rString] scanHexInt:&r];
    [[NSScanner scannerWithString:gString] scanHexInt:&g];
    [[NSScanner scannerWithString:bString] scanHexInt:&b];

    NSLog(@"%f,%f,%f", (float) r, (float) g, (float) b);
    NSLog(@"%f,%f,%f,%f",(CGFloat) r / 255.0f ,(CGFloat) g / 255.0f ,(CGFloat) b / 255.0f ,(CGFloat) 1.0f);

//    NSColor *testColor = [NSColor colorWithColorSpace:colorSpace components:SRGB];
//    NSColorSpace *colorSpace = [NSColorSpace sRGBColorSpace];
//    NSColor *bColor = [aColor colorUsingColorSpaceName:NSCalibratedRGBColorSpace];


//    CGColorRef  colorRef = CGColorCreateGenericRGB(0.75, 0.2, 0.1, 1.0);
//    NSColor *color = [NSColor colorWithColorSpace:NSCalibratedRGBColorSpace components:colorRef count:4];
//    float red = 0.5f;
//    float green = 0.2f;
//    float blue = 0.4f;
//    float alpha = 0.8f;

//    colorUsingColorSpaceName:NSCalibratedRGBColorSpace
//    NSColor *color= [NSColor colorWithDeviceRed:red green:green blue:blue alpha:alpha];


//    NSColor *color = [NSColor colorW:NSCalibratedRGBColorSpace];
//    color.gr= 0.75;
//    NSColor* color = [NSColor colorWithCalibratedRed:0.75 green:0.25 blue:0.1 alpha:1];

    NSColorSpace *colorSpace = [NSColorSpace sRGBColorSpace];


//    NSColor *color = [NSColor colorWithColorSpace:colorSpace components:SRGB count:4];


//    NSColor* color = [NSColor colorWithSRGBRed:0.1 green:0.2 blue:0.3 alpha:1.0];
//    NSColor* color = [NSColor colorWithDeviceRed:0.1 green:0.2 blue:0.3 alpha:1.0];

//    [color.colorSpace:colorSpace];
//    return color;
//    return [NSColor colorWithCalibratedRed:(CGFloat) r / 255.0f green:(CGFloat) g / 255.0f blue:(CGFloat) b / 255.0f alpha:(CGFloat) 1.0f];
//    return [NSColor colorWithDeviceRed:(CGFloat) r / 255.0f green:(CGFloat) g / 255.0f blue:(CGFloat) b / 255.0f alpha:(CGFloat) 1.0f];
//    return [NSColor colorWithCalibratedRed:r/255.0f green:g/255.0f blue:b/255.0f alpha:1.0f];
    return [NSColor colorWithSRGBRed:((float) r / 255.0f)
                               green:((float) g / 255.0f)
                                blue:((float) b / 255.0f)
                               alpha:1.0f];
//    return [NSColor colorWithSRGBRed:<#(CGFloat)red#> green:<#(CGFloat)green#> blue:<#(CGFloat)blue#> alpha:<#(CGFloat)alpha#>]
}

@end
