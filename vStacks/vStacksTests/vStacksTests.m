//
//  vStacksTests.m
//  vStacksTests
//
//  Created by Nathan Dunn on 1/6/14.
//  Copyright (c) 2014 Nathan Dunn. All rights reserved.
//

#import <XCTest/XCTest.h>

@interface vStacksTests : XCTestCase

@end

@implementation vStacksTests

- (void)setUp
{
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown
{
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

- (void)testExample
{
    NSString* testString = @"TGCAGGCCCTGGAGGAGGAGTTTTCCAGGAAGCTGCAGGAACAGGAAGTGTTCTTTAAGATGAGCGGCGAATCGGAGTGCCTTAACCCCTCCAGC" ;
    NSMutableString *testAString = [NSMutableString stringWithString:testString];
   
    NSLog(@"testString length %ld",testString.length);
    XCTAssert(testString.length==95, @"should be 94");
    XCTAssert(testAString.length==95, @"should be 94");
    
    [testAString insertString:@"test" atIndex:95];
    
    NSLog(@"testString length %ld",testString.length);
    XCTAssert(testAString.length==95+4, @"should be 98");
    
    
//    XCTFail(@"No implementation for \"%s\"", __PRETTY_FUNCTION__);
}

@end
