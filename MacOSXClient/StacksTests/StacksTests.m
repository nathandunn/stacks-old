//
//  StacksTests.m
//  StacksTests
//
//  Created by Nathan Dunn on 3/13/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StacksTests.h"
#import "StacksLoader.h"
#import "StacksDocument.h"

@implementation StacksTests

- (void)setUp
{
    [super setUp];
    
    // Set-up code here.
}

- (void)tearDown
{
    // Tear-down code here.
    
    [super tearDown];
}

- (void)testLoadLoci
{
//    STFail(@"Unit tests are not implemented yet in StacksTests");
    STAssertNil(nil,@"should obviously be nil");

    StacksLoader *stacksLoader = [[StacksLoader alloc] init];
    NSString *examplePath = @"/tmp/stacks_tut/";
    StacksDocument* stacksDocument = [stacksLoader loadLoci:examplePath];
    NSUInteger lociCount = stacksDocument.locusViews.count;
    NSLog(@"lociCount %ld",lociCount);
    STAssertTrue(462==lociCount,@"locusviews should not be empty");
}

- (void)testLoadLociAndGenotypes
{
//    STFail(@"Unit tests are not implemented yet in StacksTests");
    STAssertNil(nil,@"should obviously be nil");

    StacksLoader *stacksLoader = [[StacksLoader alloc] init];
    NSString *examplePath = @"/tmp/stacks_tut/";
    StacksDocument* stacksDocument = [stacksLoader loadLociAndGenotypes:examplePath];
    NSUInteger lociCount = stacksDocument.locusViews.count;
    NSLog(@"lociCount %ld",lociCount);
    STAssertTrue(462==lociCount,@"locusviews should not be empty");
}

@end
