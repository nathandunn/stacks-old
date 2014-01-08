//
//  vStacksTests.m
//  vStacksTests
//
//  Created by Nathan Dunn on 1/6/14.
//  Copyright (c) 2014 Nathan Dunn. All rights reserved.
//

#import <XCTest/XCTest.h>
#import "CHCSVParser.h"

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
    NSString *file = [[NSBundle bundleForClass:[self class]] pathForResource:@"test" ofType:@"tsv"];
    
    NSArray *fields = [NSArray arrayWithContentsOfTSVFile:file options:CHCSVParserOptionsRecognizesBackslashesAsEscapes];
    NSLog(@"# of rows %ld",fields.count);
    
    for(int i = 0 ; i < fields.count ; i++){
        NSArray* field = [fields objectAtIndex:i];
        if(field.count==1){
            XCTAssert([[field objectAtIndex:0] isEqualToString:@""], @"should be an empty string");
        }
        else{
            XCTAssert(field.count>1, @"should have parsed");
            XCTAssert(field.count<20, @"should have not everything");
        }
//        NSLog(@"read: %ld", field.count);
    }
    XCTAssert(fields.count>10, @"test succeeds");
    
//    NSArray *expectedFields = [self expectedFields];
//    XCTFail(@"No implementation for \"%s\"", __PRETTY_FUNCTION__);
}

@end
