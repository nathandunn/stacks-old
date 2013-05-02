//
//  StacksDocument.m
//  StacksGui3
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StacksDocument.h"

@implementation StacksDocument

// TODO: remove these in favor of NSSet loci
@synthesize locusViews;
@synthesize loci;

- (id)init
{
    self = [super init];
    if (self) {
        // Add your subclass-specific initialization here.
    }
    return self;
}

//- (void)makeWindowControllers {
//    [super makeWindowControllers];
//}


- (NSString *)windowNibName
{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"StacksDocument";
}

- (void)windowControllerDidLoadNib:(NSWindowController *)aController
{
//    NSString *aControllerName = [anIdentifier stringByAppendingString: @"ViewController"];
//    NSString *aNibName = [anIdentifier stringByAppendingString: @"View"];
//    Class aControllerClass = NSClassFromString(aControllerName);
//    [self setCurrentController: [[aControllerClass alloc] initWithNibName: aNibName bundle: [NSBundle mainBundle]]];

    [super windowControllerDidLoadNib:aController];
    // Add any code here that needs to be executed once the windowController has loaded the document's window.


}

+ (BOOL)autosavesInPlace
{
    return YES;
}

- (id)initWithLocusView:(NSMutableDictionary*)locusViews {
    if ((self = [super init])) {
        self.locusViews = locusViews;
        self.orderedLocus = [[locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
            return [obj1 integerValue] - [obj2 integerValue];
        }];

    }
    return self ;
}

- (id)initWithLoci:(NSSet*)loci {
    if ((self = [super init])) {
        self.loci = loci ;
//        self.orderedLocus = [[locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
//            return [obj1 integerValue] - [obj2 integerValue];
//        }];

    }
    return self ;
}


- (NSMutableArray *)findPopulations {

    NSMutableArray *populations = [[NSMutableArray alloc] init];

    NSString *population ;
    for(population in self.populationLookup.allValues){
        if(![populations containsObject:population]){
            [populations addObject:population];
        }
    }

    return populations;
}

- (NSManagedObjectContext *)getContext {

    NSDictionary *options = [NSDictionary dictionaryWithObjectsAndKeys:[NSNumber numberWithBool:YES],NSMigratePersistentStoresAutomaticallyOption,
                                                                       [NSNumber numberWithBool:YES],
                                                                       NSInferMappingModelAutomaticallyOption, nil];


    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
    NSString *basePath = ([paths count] > 0) ? [paths objectAtIndex:0] : nil;
    NSURL *storeUrl = [NSURL fileURLWithPath:[basePath stringByAppendingFormat:@"StacksDocument.sqlite"]];
    NSPersistentStoreCoordinator *persistentStoreCoordinator = [[NSPersistentStoreCoordinator alloc] initWithManagedObjectModel:[NSManagedObjectModel mergedModelFromBundles:nil]];
    NSError *error = nil;

    if (![persistentStoreCoordinator addPersistentStoreWithType:NSSQLiteStoreType configuration:nil URL:storeUrl options:options error:&error]) {
        NSLog(@"error loading persistent store..");
        [[NSFileManager defaultManager] removeItemAtPath:storeUrl.path error:nil];
        if (![persistentStoreCoordinator addPersistentStoreWithType:NSSQLiteStoreType configuration:nil URL:storeUrl options:options error:&error]) {
            NSLog(@"Unresolved error %@, %@", error, [error userInfo]);
            //abort();
        }
    }


    NSManagedObjectContext *context = [[NSManagedObjectContext alloc] init];
    [context setPersistentStoreCoordinator:persistentStoreCoordinator];

    return context;


}

//- (NSManagedObjectContext *)managedObjectContext {
//    static NSManagedObjectContext *moc = nil;
//
//    if (moc != nil) {
//        return moc;
//    }
//
//    NSPersistentStoreCoordinator *coordinator =
//            [[NSPersistentStoreCoordinator alloc]
//                    initWithManagedObjectModel: [self managedObjectModel]];
//
//    NSString *STORE_TYPE = NSXMLStoreType;
//    NSString *STORE_FILENAME = @"CDCLI.cdcli";
//
//    NSError *error;
//    NSURL *url = [[self applicationLogDirectory] URLByAppendingPathComponent:STORE_FILENAME];
//
//    NSPersistentStore *newStore = [coordinator addPersistentStoreWithType:STORE_TYPE
//                                                            configuration:nil URL:url options:nil
//                                                                    error:&error];
//
//    if (newStore == nil) {
//
//        NSLog(@"Store Configuration Failure\n%@",
//                ([error localizedDescription] != nil) ?
//                        [error localizedDescription] : @"Unknown Error");
//    }
//
//
//    moc = [[NSManagedObjectContext alloc] initWithConcurrencyType:NSMainQueueConcurrencyType];
//    [moc setPersistentStoreCoordinator:coordinator];
//
//    return moc;
//}

- (NSManagedObjectModel*) managedObjectModel {

    static NSManagedObjectModel *mom = nil;

    if (mom != nil) {
        return mom;
    }

    NSEntityDescription *runEntity = [[NSEntityDescription alloc] init];
    [runEntity setName:@"Run"];
    [runEntity setManagedObjectClassName:@"Run"];

    NSAttributeDescription *dateAttribute = [[NSAttributeDescription alloc] init];

    [dateAttribute setName:@"date"];
    [dateAttribute setAttributeType:NSDateAttributeType];
    [dateAttribute setOptional:NO];


    NSAttributeDescription *idAttribute = [[NSAttributeDescription alloc] init];

    [idAttribute setName:@"processID"];
    [idAttribute setAttributeType:NSInteger64AttributeType];
    [idAttribute setOptional:NO];
    [idAttribute setDefaultValue:@(-1)];

    NSExpression *lhs = [NSExpression expressionForEvaluatedObject];
    NSExpression *rhs = [NSExpression expressionForConstantValue:@0];

    NSPredicate *validationPredicate = [NSComparisonPredicate
            predicateWithLeftExpression:lhs
                        rightExpression:rhs
                               modifier:NSDirectPredicateModifier
                                   type:NSGreaterThanPredicateOperatorType
                                options:0];

    NSString *validationWarning = @"Process ID < 1";

    [idAttribute setValidationPredicates:@[validationPredicate]
                  withValidationWarnings:@[validationWarning]];

    [runEntity setProperties:@[dateAttribute, idAttribute]];

    mom = [[NSManagedObjectModel alloc] init];
    [mom setEntities:@[runEntity]];

    NSDictionary *localizationDictionary = @{
            @"Property/date/Entity/Run":@"Date",
            @"Property/processID/Entity/Run":@"Process ID",
            @"ErrorString/Process ID < 1":@"Process ID must not be less than 1"};

    [mom setLocalizationDictionary:localizationDictionary];

    return mom;
}

- (NSURL*) applicationLogDirectory {

    NSString *LOG_DIRECTORY = @"CDCLI";
    static NSURL *ald = nil;

    if (ald == nil) {

        NSFileManager *fileManager = [[NSFileManager alloc] init];
        NSError *error;
        NSURL *libraryURL = [fileManager URLForDirectory:NSLibraryDirectory inDomain:NSUserDomainMask appropriateForURL:nil create:YES error:&error];
        if (libraryURL == nil) {
            NSLog(@"Could not access Library directory\n%@", [error localizedDescription]);
        }
        else {
            ald = [libraryURL URLByAppendingPathComponent:@"Logs"];
            ald = [ald URLByAppendingPathComponent:LOG_DIRECTORY];
            NSDictionary *properties = [ald resourceValuesForKeys:@[NSURLIsDirectoryKey]
                                                            error:&error];
            if (properties == nil) {
                if (![fileManager createDirectoryAtURL:ald withIntermediateDirectories:YES attributes:nil error:&error]) {
                    NSLog(@"Could not create directory %@\n%@", [ald path], [error localizedDescription]);
                    ald = nil;
                }
            }
        }
    }
    return ald;
}


@end
