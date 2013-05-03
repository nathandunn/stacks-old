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

//- (NSManagedObjectContext *)getContext {
//    NSDictionary *options = [NSDictionary dictionaryWithObjectsAndKeys:[NSNumber numberWithBool:YES],NSMigratePersistentStoresAutomaticallyOption,
//                                                                       [NSNumber numberWithBool:YES],
//                                                                       NSInferMappingModelAutomaticallyOption, nil];
//
//
//    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
//    NSString *basePath = ([paths count] > 0) ? [paths objectAtIndex:0] : nil;
//    NSURL *storeUrl = [NSURL fileURLWithPath:[basePath stringByAppendingFormat:@"/StacksDocument.sqlite"]];
//    NSLog(@"saving to %@ from %@",basePath,storeUrl);
//    NSPersistentStoreCoordinator *persistentStoreCoordinator = [[NSPersistentStoreCoordinator alloc] initWithManagedObjectModel:[NSManagedObjectModel mergedModelFromBundles:nil]];
//    NSError *error = nil;
//
//    if (![persistentStoreCoordinator addPersistentStoreWithType:NSSQLiteStoreType configuration:nil URL:storeUrl options:options error:&error]) {
//        NSLog(@"error loading persistent store..");
//        [[NSFileManager defaultManager] removeItemAtPath:storeUrl.path error:nil];
//        if (![persistentStoreCoordinator addPersistentStoreWithType:NSSQLiteStoreType configuration:nil URL:storeUrl options:options error:&error]) {
//            NSLog(@"Unresolved error %@, %@", error, [error userInfo]);
//            //abort();
//        }
//    }
//
//
//    NSManagedObjectContext *context = [[NSManagedObjectContext alloc] init];
//    [context setPersistentStoreCoordinator:persistentStoreCoordinator];
//
//    return context;
//}


- (NSManagedObjectContext *)getContextForPath:(NSString *)path {
    return [self getContextForPath:path andName:@"StacksDocument"];
}

- (NSManagedObjectContext *)getContextForPath:(NSString *)path andName:(NSString *)name {
    NSDictionary *options = [NSDictionary dictionaryWithObjectsAndKeys:[NSNumber numberWithBool:YES],NSMigratePersistentStoresAutomaticallyOption,
                                                                       [NSNumber numberWithBool:YES],
                                                                       NSInferMappingModelAutomaticallyOption, nil];


//    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
//    NSString *basePath = ([paths count] > 0) ? [paths objectAtIndex:0] : nil;
    NSURL *storeUrl = [NSURL fileURLWithPath:[path stringByAppendingFormat:@"/%@.sqlite",name]];
    NSLog(@"saving to %@ from %@",path,storeUrl);
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
@end
