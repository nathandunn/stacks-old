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

@implementation StacksGui3Tests {
//    stacksConverter;
    StacksConverter *stacksConverter ;
}


- (void)setUp {
    [super setUp];
    stacksConverter = [[StacksConverter alloc] init];
    // Set-up code here.
}

- (void)tearDown {
    // Tear-down code here.
    stacksConverter = nil ;

    [super tearDown];
}

//- (void)testExample
//{
//    STFail(@"Unit tests are not implemented yet in StacksGui3Tests");
//}

- (void)testReadRawStacks {
//    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *examplePath = @"/tmp/stacks_tut/";
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];
    if (existsAtPath) {
//        [self loadApplication:examplePath];
        StacksDocument *stacksDocument = [stacksConverter loadLociAndGenotypes:examplePath];
        NSSet *loci = stacksDocument.loci;
        STAssertEquals( (NSUInteger) 462, loci.count, @"should match loci count");


        LocusMO *locusMO = [loci.allObjects objectAtIndex:0];
        NSLog(@"locus %@ has %ld datums", locusMO.locusId, locusMO.datums.count);

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
//            NSLog(@"locus %@ has %ld datums",locusMO.locusId,locusMO.datums.count);
//        }
    }
    else {
        NSLog(@"%@ does not exist.", examplePath);
        STFail(@"Does not exists at path!");
    }
}

//- (void)testCreatePopulatedStoreToPath {
//    NSString *examplePath = @"/tmp/stacks_tut/";
//    NSString *filePath = [examplePath stringByAppendingString:@"/StacksDocument.sqlite"];
//
//    NSFileManager *fileManager = [NSFileManager defaultManager];
//    NSError *fileError;
//    if ([fileManager fileExistsAtPath:filePath]) {
//        [fileManager removeItemAtPath:filePath error:&fileError];
//    }
//    if (fileError) {
//        STFail(@"error deleting file %@", fileError);
//    }
//    STAssertFalse([fileManager fileExistsAtPath:filePath], @"Should be false");
//
//    StacksDocument *stacksDocument = [stacksConverter loadLociAndGenotypes:examplePath];
//
//    NSLog(@"loci count %ld",stacksDocument.loci.count);
//
//    STAssertTrue(stacksDocument.loci.count > 8, @"should be at least 8 loci %ld", stacksDocument.loci.count );
//
//}

- (void)testCreatePopulatedStoreToPathAndContext {
//    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
    NSString *examplePath = @"/tmp/stacks_tut/";
    NSString *filePath = [examplePath stringByAppendingString:@"/StacksDocument.sqlite"];

    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSError *fileError;
    if ([fileManager fileExistsAtPath:filePath]) {
        [fileManager removeItemAtPath:filePath error:&fileError];
    }
    if (fileError) {
        STFail(@"error deleting file %@", fileError);
    }
    STAssertFalse([fileManager fileExistsAtPath:filePath], @"Should be false");

    StacksDocument *stacksDocument = [stacksConverter createStacksDocumentForPath:examplePath];
    stacksDocument = [stacksConverter loadDocument:stacksDocument];
    if (stacksDocument == nil) {
        STFail(@"There was an error reading in the stacks Document ");
    }
    NSLog(@"loci count %ld",stacksDocument.loci.count);
    STAssertTrue(stacksDocument.loci.count > 8, @"should be at least 8 loci %ld", stacksDocument.loci.count );
    NSError *error2;
    if (![stacksDocument.managedObjectContext save: &error2]) {
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

//- (void)testCreateEmptyStore {
//    NSError *stacksDocumentCreateError;
//    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
////    NSString *examplePath = @"~/Desktop/stacks_tut/";
//    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
//    NSString *basePath = ([paths count] > 0) ? [paths objectAtIndex:0] : nil;
//    NSManagedObjectContext *moc = [stacksDocument getContextForPath:basePath andName:@"Empty"];
//
//    NSError *error2;
//    if (![moc save:&error2]) {
//        NSLog(@"Error while saving %@", error2);
//        STFail(@"Failed to save %@", error2);
//    }
//    else {
//        NSLog(@"SUCCESS!!!");
//    }
//}


- (void)testReadPopulatedDataStore {
    NSString *examplePath = @"/tmp/stacks_tut/";
    NSString *filePath = [examplePath stringByAppendingString:@"/StacksDocument.sqlite"];

    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL existsAtPath = [fileManager fileExistsAtPath:filePath];
    if (!existsAtPath) {
        NSLog(@"file does NOT exit! %@",filePath);
    }
//    StacksDocument *newStacksDocument = [stacksConverter createStacksDocumentForPath:examplePath];
//    StacksDocument *newStacksDocument = [stacksConverter getStacksDocumentForPath:examplePath];

    NSError *stacksDocumentCreateError ;
//    StacksDocument *newStacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
//    NSURL *fileUrl = [NSURL fileURLWithPath:filePath]
    NSURL *fileUrl = [NSURL fileURLWithPath:[examplePath stringByAppendingString:@"/StacksDocument.sqlite"]];
//    StacksDocument *newStacksDocument = [stacksConverter createStacksDocumentForPath:examplePath];
    StacksDocument *newStacksDocument = [[StacksDocument alloc]
            initWithContentsOfURL:fileUrl ofType:NSSQLiteStoreType error:&stacksDocumentCreateError];
    if(stacksDocumentCreateError){
        STFail(@"failed to load . . .error %@",stacksDocumentCreateError);
    }

    NSError *error ;
    NSManagedObjectContext *moc = newStacksDocument.managedObjectContext ;
    NSEntityDescription *entityDescription = [NSEntityDescription
            entityForName:@"Locus" inManagedObjectContext:moc];
    NSFetchRequest *request = [[NSFetchRequest alloc] init];
    [request setEntity:entityDescription];

    NSArray *locusArray = [moc executeFetchRequest:request error:&error];
    for(NSUInteger  i = 0 ; i < 5 ; i++){
        LocusMO *locusMO = [locusArray objectAtIndex:i];
        NSLog(@"index %ld locus %@",i,locusMO.locusId);
        NSLog(@"has datums %ld",locusMO.datums.count);
    }
    NSLog(@"array size %ld", locusArray.count);

    STAssertTrue(locusArray.count ==462 , @"should be at least 8 loci %ld", locusArray.count );
}

@end
