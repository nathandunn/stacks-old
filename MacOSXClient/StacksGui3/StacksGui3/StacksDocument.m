//
//  StacksDocument.m
//  StacksGui3
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StacksDocument.h"
#import "DatumMO.h"
#import "LocusMO.h"
#import "PopulationMO.h"
#import "DatumRepository.h"
#import "PopulationRepository.h"
#import "LocusRepository.h"
#import "PopulationArrayController.h"

@interface StacksDocument()

@property(weak) IBOutlet NSTableView *locusTableView;
@property(weak) IBOutlet NSTableView *populationTableView;
@property(weak) IBOutlet NSTableView *stacksTableView;
@property(weak) IBOutlet NSCollectionView *datumCollectionView;
@property(weak) IBOutlet NSArrayController *datumController ;
//@property(weak) IBOutlet PopulationArrayController *populationController ;

//@property(weak) IBOutlet NSArrayController *stacksController ;

@end

@implementation StacksDocument

// TODO: remove these in favor of NSSet loci
//@synthesize locusViews;
@synthesize loci;
@synthesize populations;

// selected stuff
@synthesize selectedDatums;
@synthesize selectedDatum;

// repository
@synthesize datumRepository;
@synthesize locusRepository;
@synthesize populationRepository;

// array controller
@synthesize datumController ;
//@synthesize populationController;

- (id)init {
    self = [super init];
    if (self) {
        // Add your subclass-specific initialization here.
        datumRepository = [[DatumRepository alloc] init];
        locusRepository = [[LocusRepository alloc] init];
        populationRepository = [[PopulationRepository alloc] init];
        
//        [datumController addObserver:self forKeyPath:@"selectionIndexes" options:NSKeyValueObservingOptionInitial context:nil];
    }
    return self;
}

//- (void)makeWindowControllers {
//    [super makeWindowControllers];
//}


- (NSString *)windowNibName {
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"StacksDocument";
}

- (void)windowControllerDidLoadNib:(NSWindowController *)aController {
//    NSString *aControllerName = [anIdentifier stringByAppendingString: @"ViewController"];
//    NSString *aNibName = [anIdentifier stringByAppendingString: @"View"];
//    Class aControllerClass = NSClassFromString(aControllerName);
//    [self setCurrentController: [[aControllerClass alloc] initWithNibName: aNibName bundle: [NSBundle mainBundle]]];
    [datumController addObserver:self forKeyPath:@"selectionIndexes" options:(NSKeyValueObservingOptionNew) context:nil];
//    [populationController addObserver:self forKeyPath:@"selectionIndexes" options:(NSKeyValueObservingOptionNew) context:nil];
//    [datumController addObserver:self forKeyPath:@"selectionIndexes" options:(NSKeyValueObservingOptionInitial | NSKeyValueObservingOptionNew | NSKeyValueObservingOptionOld) context:nil];

    [super windowControllerDidLoadNib:aController];
    // Add any code here that needs to be executed once the windowController has loaded the document's window.
//    [datumController addObserver:self forKeyPath:@"selectionIndexes" options:(NSKeyValueObservingOptionNew | NSKeyValueObservingOptionOld) context:nil];
}

+ (BOOL)autosavesInPlace {
    return YES;
}

- (void)tableViewSelectionDidChange:(NSNotification *)aNotification {
//    NSString *tableName = [[aNotification object] identifier];
//    NSLog(@"table selected!! %@",tableName);

    self.selectedLocus = [self findSelectedLocus];
    self.selectedPopulation = [self findSelectedPopulation];

//    NSLog(@"selected locus: %@ and population: %@",self.selectedLocus.locusId,self.selectedPopulation.populationId);
    if(self.selectedLocus!=nil && self.selectedPopulation!=nil){
//        self.selectedDatums = [self.datumRepository getDatums:self.managedObjectContext locus:self.selectedLocus andPopulation:self.selectedPopulation];
        self.selectedDatums = [self.datumRepository getDatumsOrdered:self.managedObjectContext locus:self.selectedLocus andPopulation:self.selectedPopulation];
        self.selectedDatum = [self.selectedDatums objectAtIndex:0];
    }
    else{
        self.selectedDatums = nil ;
        self.selectedDatum = nil ;
    }
    
}

- (PopulationMO *)findSelectedPopulation {
    NSInteger selectedRow = [self.populationTableView selectedRow];
    if(selectedRow>=0){
        return [[populationRepository getAllPopulations:self.managedObjectContext] objectAtIndex:selectedRow];
    }
    return nil ;
}

- (LocusMO *)findSelectedLocus {
    NSInteger selectedRow = [self.locusTableView selectedRow];
    if(selectedRow>=0){
        return [locusRepository getLocus:self.managedObjectContext forId:selectedRow];
//        return (LocusMO*) [[self.loci allObjects] objectAtIndex: (NSUInteger) selectedRow];
    }
    return nil ;
//    NSArray *sortedKeys = [[self.stacksDocument.locusViews allKeys] sortedArrayUsingComparator:(NSComparator) ^(id obj1, id obj2) {
//        return [obj1 integerValue] - [obj2 integerValue];
//    }];
//    NSString *key = [sortedKeys objectAtIndexedSubscript:selectedRow];
//    LocusView *locusView = [self.stacksDocument.locusViews objectForKey:key];
//    return locusView;
}


//- (NSMutableArray *)findPopulations {
//
//    NSMutableArray *populations = [[NSMutableArray alloc] init];
//
//    NSString *population;
//    for (population in self.populationLookup.allValues) {
//        if (![populations containsObject:population]) {
//            [populations addObject:population];
//        }
//    }
//
//    return populations;
//}

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
    NSDictionary *options = [NSDictionary dictionaryWithObjectsAndKeys:[NSNumber numberWithBool:YES], NSMigratePersistentStoresAutomaticallyOption,
                                                                       [NSNumber numberWithBool:YES],
                                                                       NSInferMappingModelAutomaticallyOption, nil];


    NSURL *storeUrl = [NSURL fileURLWithPath:[path stringByAppendingFormat:@"/%@.sqlite", name]];
    NSLog(@"saving to %@ from %@", path, storeUrl);
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
    self.managedObjectContext = context;
    self.path = path;


    return context;
}


- (BOOL)readFromURL:(NSURL *)absoluteURL ofType:(NSString *)typeName error:(NSError **)error {
    BOOL returnType = [super readFromURL:absoluteURL ofType:typeName error:error];
//    NSLog(@"reading from the url %@", absoluteURL);
//    if(error!=nil){
//        NSLog(@"error reading %@", error);
//    }
//    NSLog(@"error reading success %ld",  returnType);

    return returnType;
}

- (void)observeValueForKeyPath:(NSString *)keyPath
                      ofObject:(id)object
                        change:(NSDictionary *)change
                       context:(void *)context {
    if ([keyPath isEqualTo:@"selectionIndexes"]) {
        if ([[datumController selectedObjects] count] > 0) {
            if ([[datumController selectedObjects] count] == 1) {
                DatumMO *datumMO = (DatumMO *) [[datumController selectedObjects] objectAtIndex:0];
                if ([datumMO.name isEqualToString:self.previousStacksName]) {
                    return;
                }
                self.selectedDatum = datumMO;
                self.previousStacksName = datumMO.name ;
            }

        }

    }
}

@end
