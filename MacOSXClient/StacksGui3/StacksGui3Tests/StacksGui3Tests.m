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
        NSURL *storeURL = [NSURL URLWithString:[filePath stringByAppendingFormat:@"storedNew.sqlite"]];
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


- (void)testCreateStore {
    NSError *stacksDocumentCreateError ;
    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
    NSManagedObjectContext *moc = stacksDocument.managedObjectContext;
    NSPersistentStoreCoordinator *psc = [moc persistentStoreCoordinator];
    NSError *error = nil;
    NSDictionary *options =
            [NSDictionary dictionaryWithObject:[NSNumber numberWithBool:1]
                                        forKey:NSReadOnlyPersistentStoreOption];

    NSString *examplePath = @"/Users/ndunn/Desktop/stacks_tut/stored.sqlite";
//    NSString *examplePath = @"~/Desktop/stacks_tut/";
//    NSURL *storeURL = [NSURL URLWithString:[examplePath stringByAppendingFormat:@"stored.sqlite"]];
    NSURL *storeURL = [NSURL fileURLWithPath:examplePath];
//    NSLog(@"store URL %@ fileUrl %@",storeURL,[NSURL fileURLWithPath:storeURL]);
    NSLog(@"store URL %@", storeURL);
//    NSPersistentStore *roStore =
    if (![psc addPersistentStoreWithType:NSSQLiteStoreType
                           configuration:nil URL:storeURL
                                 options:options error:&error]) {
        STFail(@"Failed to add persistent store to coordinator! %@",error);
    }
    NSError *error2;
    BOOL saved = [stacksDocument.managedObjectContext save:&error2];
    NSLog(@"saved %d error %@", saved, error2);

}

- (void)testCreateRawStore {

//    NSString *storePath = [[self applicationDocumentsDirectory] stringByAppendingPathComponent: @"Datafile.sqlite"];
    NSString *storePath = @"file://localhost/tmp/stacks_tut/StacksGui3.sqlite" ;
//    NSString *storePath = @"Store.sqlite" ;
    

    NSFileManager *fileManager = [NSFileManager defaultManager];
//    if (![fileManager fileExistsAtPath:storePath]) {
//        NSLog(@"does not exist so createing");
//        NSString *defaultStorePath = [[NSBundle mainBundle]  pathForResource:@"Datafile-DefaultData" ofType:@"sqlite"];
//        if (defaultStorePath) {
//            NSLog(@"default DOES exist") ;
//            [fileManager copyItemAtPath:defaultStorePath toPath:storePath error:NULL];
//        }
//        else{
//            NSLog(@"default DOES NOT exist") ;
//        }
//    }
//    else{
//        NSLog(@"file DOES exist") ;
//    }

    NSError *stacksDocumentCreateError ;
    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
//    NSManagedObjectModel *mom = [stacksDocument managedObjectModel];
    NSManagedObjectContext *moc = [stacksDocument getContext];
    NSPersistentStoreCoordinator *psc = [moc persistentStoreCoordinator];
    NSError *error = nil;
    NSDictionary *options =
            [NSDictionary dictionaryWithObject:[NSNumber numberWithBool:1]
                                        forKey:NSReadOnlyPersistentStoreOption];

    NSError *error1 = nil;

//    NSURL *storeURL = [NSURL fileURLWithPath:storePath];
//    NSURL *storeURL =[NSURL URLWithString:storePath];
    //    NSPersistentStore *roStore =
//    if (![psc addPersistentStoreWithType:NSSQLiteStoreType
//                           configuration:nil URL:storeURL
//                                 options:options error:&error1]) {
//        STFail(@"Failed to add persistent store to coordinator! %@",error1);
//    }

    NSError *error2;
    if (![moc save: &error2]) {
        NSLog(@"Error while saving %@",error2);
        STFail(@"Failed to save %@",error2);

//                ([error localizedDescription] != nil) ? [error localizedDescription] : @"Unknown Error");
//        exit(1);
    }
    else{
        NSLog(@"SUCCESS!!!") ;
    }




//    NSManagedObjectContext *moc = stacksDocument.managedObjectContext;
//    NSPersistentStoreCoordinator *psc = [moc persistentStoreCoordinator];
//    NSError *error = nil;
//    NSDictionary *options =
//            [NSDictionary dictionaryWithObject:[NSNumber numberWithBool:1]
//                                        forKey:NSReadOnlyPersistentStoreOption];
//
////    NSString *examplePath = @"/Users/ndunn/Desktop/stacks_tut/stored.sqlite";
//    NSString *examplePath = @"stored.sqlite";
//    //    NSString *examplePath = @"~/Desktop/stacks_tut/";
////    NSURL *storeURL = [NSURL URLWithString:[examplePath stringByAppendingFormat:@"stored.sqlite"]];
//    NSURL *storeURL = [NSURL fileURLWithPath:examplePath];
////    NSLog(@"store URL %@ fileUrl %@",storeURL,[NSURL fileURLWithPath:storeURL]);
//    NSLog(@"store URL %@", storeURL);
////    NSPersistentStore *roStore =
//    if (![psc addPersistentStoreWithType:NSSQLiteStoreType
//                           configuration:nil URL:storeURL
//                                 options:options error:&error]) {
//        STFail(@"Failed to add persistent store to coordinator! %@",error);
//    }
//    NSError *error2;
//    BOOL saved = [stacksDocument.managedObjectContext save:&error2];
//    NSLog(@"saved %d error %@", saved, error2);

}

@end
