//
//  StacksGui3Tests.m
//  StacksGui3Tests
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StacksGui3Tests.h"
#import "StacksConverter.h"
#import "StacksDocument.h"

@implementation StacksGui3Tests

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

//- (void)testExample
//{
//    STFail(@"Unit tests are not implemented yet in StacksGui3Tests");
//}

- (void)testReadRawStacks
{
    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *examplePath = @"/tmp/stacks_tut/";
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];
    if(existsAtPath){
//        [self loadApplication:examplePath];
        StacksDocument *stacksDocument= [stacksConverter loadLociAndGenotypes:examplePath];
        NSSet* loci = stacksDocument.loci;
        STAssertEquals( (NSUInteger) 462, loci.count, @"should match loci count");
    }
    else{
        NSLog(@"%@ does not exist.",examplePath);
        STFail(@"Does not exists at path!");
    }

}

@end
