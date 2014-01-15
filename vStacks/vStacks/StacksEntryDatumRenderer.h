//
// Created by Nathan Dunn on 1/14/14.
// Copyright (c) 2014 Nathan Dunn. All rights reserved.
//


#import <Foundation/Foundation.h>


@interface StacksEntryDatumRenderer : NSObject

@property(nonatomic, copy) NSMutableDictionary*snpLocusLookup;
@property(nonatomic, copy) NSMutableDictionary*snpDatumLookup;
@property(nonatomic, copy) NSNumberFormatter *numberFormatter;

- (NSString *)renderHtmlForData:(NSData *)data datumSnps:(NSData *)snps locusSnps:(NSData *)snps1;
@end