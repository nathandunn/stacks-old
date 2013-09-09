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
#import "DatumArrayController.h"
#import "StacksConverter.h"
#import "StacksDocumentController.h"

@interface StacksDocument()

@property(weak) IBOutlet NSTableView *locusTableView;
@property(weak) IBOutlet NSTableView *populationTableView;
@property(weak) IBOutlet NSTableView *stacksTableView;
@property(weak) IBOutlet NSCollectionView *datumCollectionView;
@property(weak) IBOutlet DatumArrayController *datumController ;
@property(weak) IBOutlet NSProgressIndicator *loadProgress;
@property(weak) IBOutlet NSPanel *progressPanel ;
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
@synthesize loadProgress;
@synthesize progressPanel;



- (id)init {
    self = [super init];
    if (self) {
        // Add your subclass-specific initialization here.
        datumRepository = [[DatumRepository alloc] init];
        locusRepository = [[LocusRepository alloc] init];
        populationRepository = [[PopulationRepository alloc] init];
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
    [datumController addObserver:self forKeyPath:@"selectionIndexes" options:(NSKeyValueObservingOptionNew) context:nil];

    [super windowControllerDidLoadNib:aController];
    // Add any code here that needs to be executed once the windowController has loaded the document's window.
    [self.stacksTableView setIntercellSpacing:NSMakeSize(0, 0)];
}

+ (BOOL)autosavesInPlace {
    return YES;
}



- (void)tableViewSelectionDidChange:(NSNotification *)aNotification {
    NSString *tableName = [[aNotification object] identifier];
    NSLog(@"table selected!! %@",tableName);


    self.selectedLocus = [self findSelectedLocus];
    self.selectedPopulation = [self findSelectedPopulation];

    if(self.selectedLocus!=nil && self.selectedPopulation!=nil){
        NSLog(@"getting selected locus %@",self.selectedLocus.locusId);
        NSLog(@"getting selected population %@",self.selectedPopulation.name);
        self.selectedDatums = [self.datumRepository getDatumsOrdered:self.managedObjectContext locus:self.selectedLocus andPopulation:self.selectedPopulation];
        if(self.selectedDatums!=nil && self.selectedDatums.count>0){
                self.selectedDatum = [self.selectedDatums objectAtIndex:0];
        }
        else{
            self.selectedDatum = nil ;
        }
    }
    else{
        self.selectedDatums = nil ;
        self.selectedDatum = nil ;
    }

}

- (PopulationMO *)findSelectedPopulation {
    NSInteger selectedRow = [self.populationTableView selectedRow];
    if(selectedRow>=0){
         return [populationRepository getPopulation:self.managedObjectContext byIndexSortedByName:selectedRow];
    }
    return nil ;
}

- (LocusMO *)findSelectedLocus {
    NSInteger selectedRow = [self.locusTableView selectedRow];
    if(selectedRow>=0){
        // id starts at 1 + row  . . . I hope this is always true
        return [locusRepository getLocus:self.managedObjectContext forId:selectedRow+1];
    }
    return nil ;
}


- (NSManagedObjectContext *)getContextForPath:(NSString *)path {
    return [self getContextForPath:path andName:@"StacksDocument"];
}

- (NSManagedObjectContext *)getContextForPath:(NSString *)path andName:(NSString *)name {
    NSDictionary *options = [NSDictionary dictionaryWithObjectsAndKeys:[NSNumber numberWithBool:YES], NSMigratePersistentStoresAutomaticallyOption,
                                                                       [NSNumber numberWithBool:YES],
                                                                       NSInferMappingModelAutomaticallyOption, nil];


    NSURL *storeUrl = [NSURL fileURLWithPath:[path stringByAppendingFormat:@"/%@.stacks", name]];
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

- (IBAction) importDocument:(id)sender{
    NSLog(@"Importing doc");
    NSOpenPanel * panel = [NSOpenPanel openPanel];
    [panel setAllowsMultipleSelection:NO];
    [panel setCanChooseDirectories:YES];
    [panel setCanChooseFiles:NO];
    [panel setFloatingPanel:YES];
    NSSize minSize ;
    minSize.height=600;
    minSize.width=500;

    [panel setMinSize:minSize];
    NSInteger result = [panel runModal];

    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
//    NSInteger result = [panel runModalForDirectory:NSHomeDirectory() file:nil types:nil];
    if(result == NSOKButton) {
        NSLog(@"ok !!");
        NSString *stacksDocumentPath = [panel.directoryURL.path stringByAppendingFormat:@"/StacksDocument.stacks"];
        BOOL fileRemoved = [[NSFileManager defaultManager] removeItemAtPath:stacksDocumentPath error:NULL];
        NSLog(@"file removed %i",fileRemoved);

//        StacksDocument *stacksDocument = [stacksConverter loadLociAndGenotypes:[panel.directoryURL.path stringByAppendingString:@"/"]];
//        [stacksConverter loadLociAndGenotypes:[panel.directoryURL.path stringByAppendingString:@"/"]];


        [self startProgressPanel:@"starting"];
        [stacksConverter loadLociAndGenotypes:[panel.directoryURL.path stringByAppendingString:@"/"] progressBar:loadProgress];
//        NSString* directoryStructure = []
//        NSString *stacksDocumentPath = [panel.directoryURL.path stringByAppendingFormat:@"/StacksDocument.stacks"];
//        NSString *[path stringByAppendingString:<#(NSString *)aString#>]
//        [stacksDocument open];
//        return [panel URLs];
        [[StacksDocumentController sharedDocumentController]openDocumentWithContentsOfURL: [NSURL fileURLWithPath:stacksDocumentPath] display:YES error:NULL];

        // if any other documents are around

        // set the file menu to enable "auto"

//        NSMenu *mainMenu = [[NSApplication sharedApplication] mainMenu];
//        NSMenu *appMenu = [[mainMenu itemAtIndex:1] submenu];
//        NSMenuItem *menuItem = [appMenu itemAtIndex:3];
//        [menuItem setEnabled:false];

//        for (NSMenuItem *item in [appMenu itemArray]) {
//            NSLog(@"%@", [item title]);
//        }
    }
    else{
            NSLog(@"NOT ok !!");
    }
//    return nil;
}



//- (void) startIndeterminateProgressPanel:(NSString *)message
//{
//    // Display a progress panel as a sheet
//    self.progressMessage = message;
//    [progressIndicator setIndeterminate: YES];
//    [progressIndicator startAnimation: self];
//    [progressCancelButton setEnabled: NO];
//    [NSApp beginSheet: progressPanel
//       modalForWindow: window
//        modalDelegate: self
//       didEndSelector: @selector(progressDidEnd: returnCode: contextInfo:)
//          contextInfo: NULL];
//}
- (void) startProgressPanel:(NSString *)message
{
    [loadProgress displayIfNeeded];
    [loadProgress setIndeterminate:false];
    [loadProgress setDisplayedWhenStopped:false];
    [loadProgress setNeedsDisplay:true];

//    [NSApp beginSheet: progressPanel
//       modalForWindow: self.windowForSheet
//        modalDelegate: self
//       didEndSelector: nil
//          contextInfo: nil];
//    [NSApp runModalForWindow: progressPanel];
//    // Dialog is up here.
//    [NSApp endSheet: progressPanel];
//    [progressPanel orderOut: self];

//    [NSApp beginSheet: progressPanel
//       modalForWindow: self.windowForSheet
//        modalDelegate: nil
//       didEndSelector: nil
//          contextInfo: nil];

    // Display a progress panel as a sheet
//    self.progressMessage = message;
//    [progressIndicator setIndeterminate: YES];
//    [progressIndicator startAnimation: self];
//    [progressCancelButton setEnabled: NO];

    // TODO: find acces to the modal window we are using
//    [NSApp beginSheet: progressPanel
//       modalForWindow: self.windowForSheet
//        modalDelegate: self
//       didEndSelector: @selector(progressDidEnd: returnCode: contextInfo:)
//          contextInfo: NULL];
}
- (void) progressDidEnd:(NSWindow *)panel returnCode:(int)returnCode contextInfo:(void *)context
{
//    xpc_connection_t connection = (xpc_connection_t)context;
//
//    if (returnCode != 0) {
//        // The cancel button was pressed.
//        NSBeep();
//    }
//
//    if (connection != NULL) {
//        // Cancel and release the anonymous connection which signals the remote
//        // service to stop, if working.
//        xpc_connection_cancel(connection);
//        xpc_release(connection);
//    }
}




- (void) stopProgressPanel
{

    [progressPanel orderOut: self];
    [NSApp endSheet: progressPanel returnCode: 0];
}

- (IBAction)cancelAction:(id)sender {

    [progressPanel orderOut: self];
    [NSApp endSheet: progressPanel returnCode: 1];
}


- (BOOL)validateUserInterfaceItem:(id <NSValidatedUserInterfaceItem>)anItem {
    NSLog(@"validating UI item in Stacks Document%@",anItem) ;
    return [super validateUserInterfaceItem:anItem];
}

- (BOOL)validateMenuItem:(NSMenuItem *)item {
    NSLog(@"validating in in Stacks Document menu item %@",item) ;
//    return [super validateUserInterfaceItem:item];
    if(item.tag==77){
        NSLog(@"should be returning true!");
        return YES;
    }
    else{
        return [super validateMenuItem:item];
    }
}

@end
