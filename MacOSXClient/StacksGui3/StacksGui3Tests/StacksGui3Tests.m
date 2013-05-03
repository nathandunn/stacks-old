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
#import "LocusMO.h"

@implementation StacksGui3Tests

- (void)setUp {
    [super setUp];

    // Set-up code here.
}

- (void)tearDown {
    // Tear-down code here.

    [super tearDown];
}

//- (void)testExample
//{
//    STFail(@"Unit tests are not implemented yet in StacksGui3Tests");
//}

- (void)testReadRawStacks {
    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *examplePath = @"/tmp/stacks_tut/";
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];
    if (existsAtPath) {
//        [self loadApplication:examplePath];
        StacksDocument *stacksDocument = [stacksConverter loadLociAndGenotypes:examplePath];
        NSSet *loci = stacksDocument.loci;
        STAssertEquals( (NSUInteger) 462, loci.count, @"should match loci count");


        LocusMO *locusMO = [loci.allObjects objectAtIndex:0];
        NSLog(@"locus %@ has %ld genotypes", locusMO.locusId, locusMO.genotypes.count);

//        NSURL *storeURL = <#URL for path to global store#>; // just same url
//        NSURL *storeURL = [NSURL URLWithString:[examplePath stringByAppendingFormat:@"stored.sqlite"] ];
//        id globalStore = [[stacksDocument.managedObjectContext persistentStoreCoordinator] persistentStoreForURL:storeURL];
        //        NSManagedObject *newEmployee = [NSEntityDescription
//                insertNewObjectForEntityForName:@"Employee"
//                         inManagedObjectContext:stacksDocument.managedObjectContext];
//        LocusMO *locusMO = [stacksDocument.loci.allObjects objectAtIndex:0]
//        [stacksDocument.managedObjectContext assignObject:locusMO toPersistentStore:globalStore];

//        StacksDocument *stacksDocument = [[StacksDocument alloc] init];
        NSManagedObjectContext *moc = stacksDocument.managedObjectContext;
        NSPersistentStoreCoordinator *psc = [moc persistentStoreCoordinator];
        NSDictionary *options =
                [NSDictionary dictionaryWithObject:[NSNumber numberWithBool:1]
                                            forKey:NSReadOnlyPersistentStoreOption];

        NSString *filePath = [NSString stringWithFormat:@"file://%@", examplePath];
        NSURL *storeURL = [NSURL URLWithString:[filePath stringByAppendingFormat:@"StacksDocument.sqlite"]];
//    NSLog(@"store URL %@ fileUrl %@",storeURL,[NSURL fileURLWithPath:storeURL]);
        NSLog(@"store URL %@", storeURL);

        NSError *error1 = nil;
        NSPersistentStore *roStore =
                [psc addPersistentStoreWithType:NSSQLiteStoreType
                                  configuration:nil URL:storeURL
                                        options:options error:&error1];

        NSError *error;
        BOOL saved = [stacksDocument.managedObjectContext save:&error];
        NSLog(@"saved %d error %@", saved, error);

//        for(LocusMO *locusMO in loci.allObjects){
//            NSLog(@"locus %@ has %ld genotypes",locusMO.locusId,locusMO.genotypes.count);
//        }
    }
    else {
        NSLog(@"%@ does not exist.", examplePath);
        STFail(@"Does not exists at path!");
    }
}

- (void)testCreatePopulatedStoreToPath {
    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
    NSString *examplePath = @"/tmp/stacks_tut/";
    StacksDocument *stacksDocument = [stacksConverter loadLociAndGenotypes:examplePath];
    NSManagedObjectContext *moc = [stacksDocument getContextForPath:examplePath];

    // NOT NECESSARY, but still here
    NSError *error2;
    if (![moc save: &error2]) {
        NSLog(@"Error while saving %@",error2);
        STFail(@"Failed to save %@",error2);
    }
    else{
        NSLog(@"SUCCESS!!!") ;
    }

}


- (void)testCreatePopulatedStoreToPathAndContext {
    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
    NSString *examplePath = @"/tmp/stacks_tut/";


    NSError *stacksDocumentCreateError ;
    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
    NSManagedObjectContext *moc = [stacksDocument getContextForPath:examplePath];
    stacksDocument.managedObjectContext = moc ;
    stacksDocument.path = examplePath ;
    if(stacksDocumentCreateError){
        NSLog(@"error creating stacks document %@",stacksDocumentCreateError);
        STFail(@"Failed to to create a stacks document %@",stacksDocumentCreateError);
    }

    stacksDocument = [stacksConverter loadDocument:stacksDocument];

    // NOT NECESSARY, but still here
    NSError *error2;
    if (![moc save: &error2]) {
        NSLog(@"Error while saving %@",error2);
        STFail(@"Failed to save %@",error2);
    }
    else{
        NSLog(@"SUCCESS!!!") ;
    }

}

//- (void)testCreateRawStoreWithData {
//    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
//    NSString *examplePath = @"/tmp/stacks_tut/";
//    StacksDocument *stacksDocument = [stacksConverter loadLociAndGenotypes:examplePath];
//    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
//    NSString *basePath = ([paths count] > 0) ? [paths objectAtIndex:0] : nil;
//    NSManagedObjectContext *moc = [stacksDocument getContextForPath:basePath];
//
//    NSError *error2;
//    if (![moc save: &error2]) {
//        NSLog(@"Error while saving %@",error2);
//        STFail(@"Failed to save %@",error2);
//    }
//    else{
//        NSLog(@"SUCCESS!!!") ;
//    }
//}

- (void)testCreateEmptyStore {
    NSError *stacksDocumentCreateError ;
    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
//    NSString *examplePath = @"~/Desktop/stacks_tut/";
    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
    NSString *basePath = ([paths count] > 0) ? [paths objectAtIndex:0] : nil;
    NSManagedObjectContext *moc = [stacksDocument getContextForPath:basePath andName:@"Empty"];

    NSError *error2;
    if (![moc save: &error2]) {
        NSLog(@"Error while saving %@",error2);
        STFail(@"Failed to save %@",error2);
    }
    else{
        NSLog(@"SUCCESS!!!") ;
    }
}

@end
