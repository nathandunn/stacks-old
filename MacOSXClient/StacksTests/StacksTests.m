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
#import "LocusView.h"

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
//    for( NSString *aKey in [stacksDocument.locusViews allKeys] ){
////        NSLog(@"key %@",aKey);
//        NSLog(@"key %@",aKey);
//    }

    for( LocusView *locus in [stacksDocument.locusViews allValues] ){
//        NSLog(@"key %@",aKey);
        NSLog(@"genotpye count %ld",locus.genotypes.count);
    }

    LocusView *locusView = [stacksDocument.locusViews objectForKey:[NSString stringWithFormat:@"%d",1]];
    STAssertNotNil(locusView, @"should not be nil");

    NSMutableDictionary *genotypes = locusView.genotypes;
    NSLog(@"number of genotypes: %d",genotypes.count);
//    for( NSString *aKey in [genotypes allKeys] ){
//        NSLog(@"key %@",aKey);
//    }

    
}

@end
