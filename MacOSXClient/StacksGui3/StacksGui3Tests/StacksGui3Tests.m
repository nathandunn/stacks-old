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
#import "PopulationMO.h"
#import "DatumMO.h"
#import "DatumRepository.h"
#import "LocusRepository.h"
#import "PopulationRepository.h"
#import "SampleMO.h"
#import "SampleRepository.h"

@implementation StacksGui3Tests {
//    stacksConverter;
    StacksConverter *stacksConverter;
    DatumRepository *datumRepository;
    LocusRepository *locusRepository;
    PopulationRepository *populationRepository;
    SampleRepository *sampleRepository;
}


- (void)setUp {
    [super setUp];
    stacksConverter = [[StacksConverter alloc] init];

    datumRepository = [[DatumRepository alloc] init];
    locusRepository = [[LocusRepository alloc] init];
    populationRepository = [[PopulationRepository alloc] init];
    sampleRepository = [[SampleRepository alloc] init];
}

- (void)tearDown {
    // Tear-down code here.
    stacksConverter = nil ;

    datumRepository = nil ;
    locusRepository = nil ;
    populationRepository = nil ;

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
        StacksDocument *stacksDocument = [stacksConverter loadLociAndGenotypes:examplePath progressWindow:nil ];
        NSSet *loci = stacksDocument.loci;
        STAssertEquals( (NSUInteger) 462, loci.count, @"should match loci count");


        LocusMO *locusMO = [loci.allObjects objectAtIndex:0];
        NSLog(@"locus %@ has %ld datums", locusMO.locusId, locusMO.datums.count);

        NSManagedObjectContext *moc = stacksDocument.managedObjectContext;
        NSPersistentStoreCoordinator *psc = [moc persistentStoreCoordinator];
        NSDictionary *options =
                [NSDictionary dictionaryWithObject:[NSNumber numberWithBool:1]
                                            forKey:NSReadOnlyPersistentStoreOption];

        NSString *filePath = [NSString stringWithFormat:@"file://%@", examplePath];
        NSURL *storeURL = [NSURL URLWithString:[filePath stringByAppendingFormat:@"StacksDocument.stacks"]];
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
    NSString *filePath = [examplePath stringByAppendingString:@"/StacksDocument.stacks"];

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
    NSManagedObjectContext *moc = stacksDocument.managedObjectContext;
    stacksDocument = [stacksConverter loadDocument:stacksDocument progressWindow:0];
    if (stacksDocument == nil) {
        STFail(@"There was an error reading in the stacks Document ");
    }
    NSLog(@"loci count %ld", stacksDocument.loci.count);
    STAssertTrue(stacksDocument.loci.count > 8, @"should be at least 8 loci %ld", stacksDocument.loci.count );
    NSError *error2;
    if (![moc save:&error2]) {
        NSLog(@"Error while saving %@", error2);
        STFail(@"Failed to save %@", error2);
    }
    else {
        NSLog(@"SUCCESS!!!");
    }

    if (error2 != nil) {
        STFail(@"Failed to save %@", error2);
    }

    NSArray *datumArray = [datumRepository getAllDatum:moc];
    NSLog(@"number of datum %ld", datumArray.count);
    NSRange nsRange;
    nsRange.length = 10;
    nsRange.location = 0;
    for (DatumMO *datumMO in [datumArray subarrayWithRange:nsRange]) {
//        NSLog(@"# of entries per stack %ld for sample %@ and loci %@", stackMO.stackEntries.count, stackMO.datum.sample.name, stackMO.datum.locus.locusId);
        STAssertTrue(datumMO.stackEntries.count > 5, @"should have atleast 5 %ld", datumMO.stackEntries.count);
//        NSLog(@"processing snps %@ ",datumMO.snps) ;
        if (datumMO.snps != nil && datumMO.snps.count > 0) {
            NSLog(@"snps count on datum %ld for sample %@ and locus %@", datumMO.snps.count, datumMO.sample.name, datumMO.locus.locusId);
        }
        else {
            NSLog(@"has no snps for sample %@ and locus %@", datumMO.sample.name, datumMO.locus.locusId);
        }
//        STAssertTrue(datumMO.snps.count>0, @"Should have snps on Datum, %ld", datumMO.snps.count);
    }

//    NSLog(@"datums found %ld",datums.count) ;
//    for(int i = 0 ; i < 5 ; i++){
//        DatumMO *datumMO = [datums objectAtIndex:i] ;
//        if(datumMO.snps != nil){
//            NSLog(@"snps count on datum %@ for sample %@ and locus %@",datumMO.snps.count,datumMO.sample.name,datumMO.locus.locusId);
//        }
//        else{
//            NSLog(@"has no snps for sample %@ and locus %@",datumMO.sample.name,datumMO.locus.locusId);
//        }
////        STAssertTrue(datumMO.snps.count>0, @"Should have snps on Datum, %ld", datumMO.snps.count);
//    }

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
    NSString *filePath = [examplePath stringByAppendingString:@"/StacksDocument.stacks"];

    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL existsAtPath = [fileManager fileExistsAtPath:filePath];
    if (!existsAtPath) {
        NSLog(@"file does NOT exit! %@", filePath);
    }
//    StacksDocument *newStacksDocument = [stacksConverter createStacksDocumentForPath:examplePath];
//    StacksDocument *newStacksDocument = [stacksConverter getStacksDocumentForPath:examplePath];

    NSError *stacksDocumentCreateError;
//    StacksDocument *newStacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
//    NSURL *fileUrl = [NSURL fileURLWithPath:filePath]
    NSURL *fileUrl = [NSURL fileURLWithPath:[examplePath stringByAppendingString:@"/StacksDocument.stacks"]];
//    StacksDocument *newStacksDocument = [stacksConverter createStacksDocumentForPath:examplePath];
    StacksDocument *newStacksDocument = [[StacksDocument alloc]
            initWithContentsOfURL:fileUrl ofType:NSSQLiteStoreType error:&stacksDocumentCreateError];
    if (stacksDocumentCreateError) {
        STFail(@"failed to load . . .error %@", stacksDocumentCreateError);
    }

    NSError *error;
    NSManagedObjectContext *moc = newStacksDocument.managedObjectContext;
    NSEntityDescription *entityDescription = [NSEntityDescription
            entityForName:@"Locus" inManagedObjectContext:moc];
    NSFetchRequest *request = [[NSFetchRequest alloc] init];
    [request setEntity:entityDescription];

    NSArray *locusArray = [moc executeFetchRequest:request error:&error];
    int locusSnpCount = 0;
    int locusAlleleCount = 0;
    for (NSUInteger i = 0; i < 5; i++) {
        LocusMO *locusMO = [locusArray objectAtIndex:i];
        NSLog(@"index %ld locus %@", i, locusMO.locusId);
        NSLog(@"has datums %ld", locusMO.datums.count);
        locusAlleleCount += locusMO.alleles.count;
        locusSnpCount += locusMO.snps.count;
    }
    NSLog(@"number of loci %ld", locusArray.count);

    STAssertTrue(locusArray.count == 462 || locusArray.count == 11, @"should be either 11 or 462 loci %ld", locusArray.count );
    STAssertTrue(locusSnpCount > 0, @"should be atleast one snp ", locusSnpCount );
    STAssertTrue(locusAlleleCount > 0, @"should be atleast one allele ", locusAlleleCount);


    NSEntityDescription *entityDescription2 = [NSEntityDescription
            entityForName:@"Population" inManagedObjectContext:moc];
    NSFetchRequest *request2 = [[NSFetchRequest alloc] init];
    [request2 setEntity:entityDescription2];

    NSArray *populationArray = [moc executeFetchRequest:request2 error:&error];
    NSLog(@"!!!population size %ld", populationArray.count);

    STAssertTrue(populationArray.count == 3, @"should be a population of 3: %ld", populationArray.count );

    PopulationMO *populationMO = [populationArray objectAtIndex:0];
    NSLog(@"index population %@", populationMO.name);
    NSLog(@"has samples %ld", populationMO.samples.count);


    NSArray *datumArray = [datumRepository getAllDatum:moc];
    NSRange nsRange;
    nsRange.length = 10;
    nsRange.location = 0;
    for (DatumMO *mo in [datumArray subarrayWithRange:nsRange]) {
        STAssertTrue(mo.stackEntries.count > 5, @"should have atleast 5 %ld", mo.stackEntries.count);
    }

    NSLog(@"num datum: %ld", datumArray.count);
    DatumMO *datumMO = [datumArray objectAtIndex:0];
    STAssertNotNil(datumMO.locus, @"should have a valid locus ");
    STAssertNotNil(datumMO.sample, @"should have a valid sample");
    STAssertNotNil(datumMO.name, @"should have a valid name ? ");
    NSSet *stackEntries = datumMO.stackEntries;
    STAssertTrue(stackEntries.count > 5, @"should have at least 5 entries %ld", stackEntries.count);
    STAssertTrue(stackEntries.count < 1000, @"but less than 1000 entries %ld", stackEntries.count);


    STAssertTrue(datumMO.depths.count > 0, @"should have at least one depth");

}

- (void)testReadWithRepositoryMethods {
    NSString *examplePath = @"/tmp/stacks_tut/";

    NSError *stacksDocumentCreateError;
    NSURL *fileUrl = [NSURL fileURLWithPath:[examplePath stringByAppendingString:@"/StacksDocument.stacks"]];
    StacksDocument *newStacksDocument = [[StacksDocument alloc] initWithContentsOfURL:fileUrl ofType:NSSQLiteStoreType error:&stacksDocumentCreateError];
    if (stacksDocumentCreateError) {
        STFail(@"failed to load . . .error %@", stacksDocumentCreateError);
    }

    NSManagedObjectContext *managedObjectContext = newStacksDocument.managedObjectContext;


    LocusMO *locusMO = [[locusRepository getAllLoci:managedObjectContext] objectAtIndex:0];
    PopulationMO *populationMO = [[populationRepository getAllPopulations:managedObjectContext] objectAtIndex:0];

    NSArray *datumWithPopulation = [datumRepository getDatums:managedObjectContext locus:locusMO andPopulation:populationMO];
    STAssertTrue(datumWithPopulation.count > 0, @"must have a population count greater than 0") ;
    NSArray *datums = [datumRepository getAllDatum:managedObjectContext];

    NSLog(@"datums found %ld", datums.count);
    int snpsCount = 0;
    int allelesCount = 0;
    for (int i = 0; i < 5; i++) {
        DatumMO *datumMO = [datums objectAtIndex:i];
        if (datumMO.snps != nil) {
            NSLog(@"snps count on datum %ld for sample %@ and locus %@", datumMO.snps.count, datumMO.sample.name, datumMO.locus.locusId);
            snpsCount += datumMO.snps.count;
        }
        else {
            NSLog(@"has no snps for sample %@ and locus %@", datumMO.sample.name, datumMO.locus.locusId);
        }

        if (datumMO.alleles != nil) {
            allelesCount += datumMO.alleles.count;
            NSLog(@"alleles count on datum %ld for sample %@ and locus %@", datumMO.alleles.count, datumMO.sample.name, datumMO.locus.locusId);
        }
        else {
            NSLog(@"has no alleles for sample %@ and locus %@", datumMO.sample.name, datumMO.locus.locusId);
        }
    }
    STAssertTrue(snpsCount > 0, @"Should have snps on Datum, %ld", snpsCount);
    STAssertTrue(allelesCount > 0, @"Should have alleles on Datum, %ld", allelesCount);
}


- (void)testCreateLargeStore {
//    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
    NSString *examplePath = @"/Users/NathanDunn/Desktop/stacks_large/";
    NSString *filePath = [examplePath stringByAppendingString:@"/StacksDocument.stacks"];

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
    NSManagedObjectContext *moc = stacksDocument.managedObjectContext;
    stacksDocument = [stacksConverter loadDocument:stacksDocument progressWindow:0];
    if (stacksDocument == nil) {
        STFail(@"There was an error reading in the stacks Document ");
    }
    NSLog(@"loci count %ld", stacksDocument.loci.count);
    STAssertTrue(stacksDocument.loci.count > 8, @"should be at least 8 loci %ld", stacksDocument.loci.count );
    NSError *error2;
    if (![moc save:&error2]) {
        NSLog(@"Error while saving %@", error2);
        STFail(@"Failed to save %@", error2);
    }
    else {
        NSLog(@"SUCCESS!!!");
    }

    NSArray *datumArray = [datumRepository getAllDatum:moc];
    NSRange nsRange;
    nsRange.length = 10;
    nsRange.location = 0;
    for (DatumMO *datumMO in [datumArray subarrayWithRange:nsRange]) {
//        NSLog(@"# of entries per stack %ld for sample %@ and loci %@", stackMO.stackEntries.count, stackMO.datum.sample.name, stackMO.datum.locus.locusId);
        STAssertTrue(datumMO.stackEntries.count > 5, @"should have atleast 5 %ld", datumMO.stackEntries.count);
    }

}

- (void)testReadLargeDataStore {
    NSString *examplePath = @"/Users/NathanDunn/Desktop/stacks_large/";
    NSString *filePath = [examplePath stringByAppendingString:@"/StacksDocument.stacks"];

    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL existsAtPath = [fileManager fileExistsAtPath:filePath];
    if (!existsAtPath) {
        NSLog(@"file does NOT exit! %@", filePath);
    }
//    StacksDocument *newStacksDocument = [stacksConverter createStacksDocumentForPath:examplePath];
//    StacksDocument *newStacksDocument = [stacksConverter getStacksDocumentForPath:examplePath];

    NSError *stacksDocumentCreateError;
//    StacksDocument *newStacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
//    NSURL *fileUrl = [NSURL fileURLWithPath:filePath]
    NSURL *fileUrl = [NSURL fileURLWithPath:[examplePath stringByAppendingString:@"/StacksDocument.stacks"]];
//    StacksDocument *newStacksDocument = [stacksConverter createStacksDocumentForPath:examplePath];
    StacksDocument *newStacksDocument = [[StacksDocument alloc]
            initWithContentsOfURL:fileUrl ofType:NSSQLiteStoreType error:&stacksDocumentCreateError];
    if (stacksDocumentCreateError) {
        STFail(@"failed to load . . .error %@", stacksDocumentCreateError);
    }

    NSError *error;
    NSManagedObjectContext *moc = newStacksDocument.managedObjectContext;
    NSEntityDescription *entityDescription = [NSEntityDescription
            entityForName:@"Locus" inManagedObjectContext:moc];
    NSFetchRequest *request = [[NSFetchRequest alloc] init];
    [request setEntity:entityDescription];

    NSArray *locusArray = [moc executeFetchRequest:request error:&error];
    for (NSUInteger i = 0; i < 5; i++) {
        LocusMO *locusMO = [locusArray objectAtIndex:i];
        NSLog(@"index %ld locus %@", i, locusMO.locusId);
        NSLog(@"has datums %ld", locusMO.datums.count);
    }
    NSLog(@"number of loci %ld", locusArray.count);

//    STAssertTrue(locusArray.count == 462 || locusArray.count==11, @"should be either 11 or 462 loci %ld", locusArray.count );


    NSEntityDescription *entityDescription2 = [NSEntityDescription
            entityForName:@"Population" inManagedObjectContext:moc];
    NSFetchRequest *request2 = [[NSFetchRequest alloc] init];
    [request2 setEntity:entityDescription2];

    NSArray *populationArray = [moc executeFetchRequest:request2 error:&error];
    NSLog(@"!!!population size %ld", populationArray.count);

//    STAssertTrue(populationArray.count == 3, @"should be a population of 3: %ld", populationArray.count );

    PopulationMO *populationMO = [populationArray objectAtIndex:0];
    NSLog(@"index population %@", populationMO.name);
    NSLog(@"has samples %ld", populationMO.samples.count);


    NSEntityDescription *datumEntityDescription = [NSEntityDescription entityForName:@"Datum" inManagedObjectContext:moc];
    NSFetchRequest *datumRequest = [[NSFetchRequest alloc] init];
    [datumRequest setEntity:datumEntityDescription];
    NSError *datumFetchError;
    NSArray *datumArray = [moc executeFetchRequest:datumRequest error:&datumFetchError];

    NSLog(@"num datum: %ld", datumArray.count);
    DatumMO *datumMO = [datumArray objectAtIndex:0];
    NSSet *stackEntries = datumMO.stackEntries;
    STAssertTrue(stackEntries.count > 5, @"should have at least 5 entries %ld", stackEntries.count);
    STAssertTrue(stackEntries.count < 100, @"but less than 100 entries %ld", stackEntries.count);
    STAssertTrue(datumMO.depths.count > 0, @"should have at least one depth");


}


- (void)testGetFinalPath{
   NSString* path = @"/Users/NathanDunn/stacks_tut2";
    NSString* lastComponent = path.lastPathComponent;
    NSLog(@"evaluating %@",lastComponent);
    STAssertEqualObjects(@"stacks_tut2", lastComponent, @"Should be equals");
    NSLog(@"passed!!! %@",lastComponent);
    //    STAssertTrue([path isEqualToString:lastComponent], @"should be equals" );
}


@end