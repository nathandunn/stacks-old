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
//#import "PopulationArrayController.h"
#import "DatumArrayController.h"
#import "GZIP.h"
#import "StackEntryDatumMO.h"
#import "StackEntryRenderer.h"
#import "StacksEntryDatumRenderer.h"
//#import "StacksConverter.h"
//#import "StacksDocumentController.h"
#import <WebKit/WebKit.h>

@interface StacksDocument ()

@property BOOL editingPopulation;
@property NSInteger previousSelectedItem;

@property(weak) IBOutlet NSTableView *locusTableView;
//@property(weak) IBOutlet NSTableView *stacksTableView;
@property(weak) IBOutlet NSCollectionView *datumCollectionView;
@property(weak) IBOutlet DatumArrayController *datumController;
@property(weak) IBOutlet NSTextField *filteredLoci;
@property(weak) IBOutlet NSTextField *totalLoci;
@property(weak) IBOutlet NSPopUpButton *populationSelector;
@property(weak) IBOutlet NSButton *editPopulationButton;
@property(weak) IBOutlet NSTextField *populationNameField;
@property(weak) IBOutlet NSTextField *maxLocusTextField;
@property(weak) IBOutlet NSPopUpButton *maxSnpPopupButton;
@property(weak) IBOutlet WebView *stacksWebView;


@end

@implementation StacksDocument

// TODO: remove these in favor of NSSet loci
//@synthesize locusViews;
@synthesize filteredLoci;
@synthesize loci;
@synthesize populations;

// selected stuff
@synthesize selectedDatums;
@synthesize selectedDatum;

// repository
//@synthesize datumRepository;
//@synthesize [LocusRepository sharedInstance];
//@synthesize populationRepository;

// array controller
@synthesize datumController;
@synthesize totalLoci;
@synthesize populationSelector;
@synthesize populationNameField;
@synthesize editPopulationButton;


@synthesize editingPopulation;
@synthesize previousSelectedItem;
@synthesize lociLocations;
@synthesize maxLocation;
@synthesize maxLocusTextField;
@synthesize snpFilterValues;
@synthesize maxSnpPopupButton;
@synthesize stacksWebView;
//@synthesize populationController;
//@synthesize loadProgress;
//@synthesize progressPanel;



- (id)init {
    self = [super init];
    if (self) {
        // Add your subclass-specific initialization here.
//        datumRepository = [[DatumRepository alloc] init];
//        [LocusRepository sharedInstance] = [[LocusRepository alloc] init];
//        populationRepository = [[PopulationRepository alloc] init];
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
//    [self.stacksTableView setIntercellSpacing:NSMakeSize(0, 0)];
//    [self.stacksTableView setEnabled:true];


    NSInteger lociCount = [[LocusRepository sharedInstance] getAllLoci:self.managedObjectContext].count;
    NSString *newString = [NSString stringWithFormat:@"%ld", lociCount];
    [totalLoci setStringValue:newString];

    editingPopulation = false ;

    double maxLocationVariable = [[LocusRepository sharedInstance] getMaxLocation:self.managedObjectContext] / 1000000;
    maxLocusTextField.stringValue = [NSString stringWithFormat:@"%1.2f", maxLocationVariable];

    [maxSnpPopupButton selectItemAtIndex:[self getSnpFilterValues].count - 1];


    NSNumberFormatter *numberFormatter = [[NSNumberFormatter alloc] init];
    [numberFormatter setNumberStyle:NSNumberFormatterNoStyle];
//    [filteredLoci setFormatter:numberFormatter];
    [[filteredLoci cell] setFormatter:numberFormatter];

    [self updateStacksView];
}

- (void)updateStacksView {

//    DatumMO *datumMO = [[datumRepository getAllDatum:[self managedObjectContext]] objectAtIndex:0];
    if (self.selectedDatum != nil) {
//        NSLog(@"loading data %@ with url %@",self.selectedDatum.stackData, [[NSBundle mainBundle] bundleURL]);
//        [[stacksWebView mainFrame] loadHTMLString:self.selectedDatum.stackData baseURL:[[NSBundle mainBundle] bundleURL]];
//        [[stacksWebView mainFrame] loadHTMLString:self.selectedDatum.stackData baseURL:[[NSBundle mainBundle] bundleURL]];
//        [[stacksWebView mainFrame] loadData:[self.selectedDatum.stackData gunzippedData] MIMEType:@"text/html" textEncodingName:@"UTF8" baseURL:[[NSBundle mainBundle] bundleURL]];

        StackEntryDatumMO *stackEntryDatumMO = [[DatumRepository sharedInstance] getStackEntryDatum:self.managedObjectContext datum:self.selectedDatum];
        if (stackEntryDatumMO != nil && stackEntryDatumMO.stackData!=nil) {

            NSData *jsonData = [stackEntryDatumMO.stackData gunzippedData];
            NSError *error;
//            NSDictionary *stackEntryData = [NSJSONSerialization JSONObjectWithData:jsonData options:kNilOptions error:&error];

            StacksEntryDatumRenderer *stacksEntryDatumRenderer = [[StacksEntryDatumRenderer alloc] init];
            NSString *html = [stacksEntryDatumRenderer renderHtmlForData:jsonData datumSnps:self.selectedDatum.snpData locusSnps:self.selectedLocus.snpData];

//            [[stacksWebView mainFrame] loadData:jsonData MIMEType:@"text/html" textEncodingName:@"UTF8" baseURL:[[NSBundle mainBundle] bundleURL]];
            [[stacksWebView mainFrame] loadHTMLString:html baseURL:[[NSBundle mainBundle] bundleURL]];

            [stacksWebView setHidden:NO];
        }
        else {

            [stacksWebView setHidden:YES];
        }

    }
    else {
        [stacksWebView setHidden:YES];
    }

//    NSURL *url = [NSURL URLWithString:@"http://www.apple.com"];
//    NSURLRequest *urlRequest = [NSURLRequest requestWithURL:url];
////    [[[self webView] mainFrame] loadRequest:urlRequest];
//    [[stacksWebView mainFrame] loadRequest:urlRequest];
}

+ (BOOL)autosavesInPlace {
    return YES;
}

- (IBAction)updateSelections:(id)sender {
    [self tableViewSelectionDidChange:nil];
}

- (IBAction)togglePopulationEdit:(id)sender {
    NSLog(@"editing %d", editingPopulation);
    if (editingPopulation) {
        NSLog(@"setting to edit");
        editPopulationButton.title = @"Edit";

        [populationSelector setHidden:false];
        [populationNameField setHidden:true];

        [populationSelector selectItemAtIndex:previousSelectedItem];
        NSLog(@"selected item index %ld", populationSelector.indexOfSelectedItem);
    }
    else {
        NSLog(@"setting to DONE");
        previousSelectedItem = populationSelector.indexOfSelectedItem;
        editPopulationButton.title = @"Done";
        [populationSelector setHidden:true];
        [populationNameField setHidden:false];
    }

    editingPopulation = !editingPopulation;
}


- (void)tableViewSelectionDidChange:(NSNotification *)aNotification {
//    NSString *tableName = [[aNotification object] identifier];
//    NSLog(@"table selected!! %@",tableName);


    self.selectedLocus = [self findSelectedLocus];
    self.selectedPopulation = [self findSelectedPopulation];

    if (self.selectedLocus != nil && self.selectedPopulation != nil) {
//        NSLog(@"getting selected locus %@", self.selectedLocus.locusId);
//        NSLog(@"getting selected population %@", self.selectedPopulation.name);
        self.selectedDatums = [[DatumRepository sharedInstance] getDatumsOrdered:self.managedObjectContext locus:self.selectedLocus.locusId andPopulation:self.selectedPopulation];
        if (self.selectedDatums != nil && self.selectedDatums.count > 0) {
            self.selectedDatum = [self.selectedDatums objectAtIndex:0];
        }
        else {
            self.selectedDatum = nil ;
        }
    }
    else {
        self.selectedDatums = nil ;
        self.selectedDatum = nil ;
    }
    [self updateStacksView];

}

- (PopulationMO *)findSelectedPopulation {
    NSInteger selectedRow = [self.populationSelector indexOfSelectedItem];
    if (selectedRow >= 0) {
        return [[PopulationRepository sharedInstance] getPopulation:self.managedObjectContext byIndexSortedByName:selectedRow];
    }
    return nil;
}

- (LocusMO *)findSelectedLocus {
    NSInteger selectedRow = [self.locusTableView selectedRow];
    if (selectedRow >= 0) {
        // id starts at 1 + row  . . . I hope this is always true
        return [[LocusRepository sharedInstance] getLocus:self.managedObjectContext forId:selectedRow + 1];
    }
    return nil;
}


- (NSArray *)generateLociLocations {
    NSArray *locusArray = [[LocusRepository sharedInstance] getLociLocations:self.managedObjectContext];
//    NSLog(@"size of array %ld",locusArray.count);

//    if(lociLocations==nil){
//        NSMutableSet* set = [[NSMutableSet alloc] init];
//        for(LocusMO *locus in locusArray){
//            [set addObject:locus.chromosome];
//        }

//        lociLocations = [[NSSet alloc] initWithSet:set];
//        lociLocations = set ;
//    }

//    return lociLocations;
//    return set ;

    return locusArray;

//    NSMutableArray *bots = [[NSMutableArray alloc] init];
//    [bots addObject:@"dogs"];
//    [bots addObject:@"cats"];
//
//    return bots ;
}

- (NSUInteger)getMaxLocation {
    return [[LocusRepository sharedInstance] getMaxLocation:self.managedObjectContext];
}

- (BOOL)noLociLocations {
    NSUInteger locusCount = [[LocusRepository sharedInstance] getLociWithChromsomes:self.managedObjectContext].count;
//    NSLog(@"NO LOCI LOCATIONS count %ld",locusCount);
    return locusCount == 0;
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
//                NSLog(@"current name: %@ vs previous: %@",datumMO.name,self.previousStacksName);
                if (self.previousStacksName != nil && [datumMO.name isEqualToString:self.previousStacksName]) {
                    return;
                }
                self.selectedDatum = datumMO;
                self.previousStacksName = datumMO.name;

            }
        }
    }
    [self updateStacksView];
}


- (BOOL)validateUserInterfaceItem:(id <NSValidatedUserInterfaceItem>)anItem {
    NSLog(@"validating UI item in Stacks Document%@", anItem);
    return [super validateUserInterfaceItem:anItem];
}

- (BOOL)validateMenuItem:(NSMenuItem *)item {
    NSLog(@"validating in in Stacks Document menu item %@", item);
//    return [super validateUserInterfaceItem:item];
    if (item.tag == 77) {
        NSLog(@"should be returning true!");
        return YES;
    }
    else {
        return [super validateMenuItem:item];
    }
}

- (void)showHelp:(id)sender {
    NSLog(@"coalling stacks doc help");
//    NSString *help = [[NSBundle mainBundle] pathForResource:@"Some Help" ofType:@"html"];
//    NSString *help = [[NSBundle mainBundle] pathForResource:@"Some Help" ofType:@"html"];
//    NSURL *url = [NSURL :@"http://creskolab.uoregon.edu/stacks/index.html"];
    NSString *url = @"http://creskolab.uoregon.edu/stacks";
//    NSLog(@"%@", url);
//    [[NSApplication sharedApplication]
//    NSString* url = @"NSLog(@\"%@\", [NSURL URLWithString:[self.storyLink description]])";
    [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:url]];
}

- (NSArray *)getSnpFilterValues {
    if (snpFilterValues == nil) {
        snpFilterValues = [NSMutableArray array];
        for (int i = 0; i < 15; i++) {
            [snpFilterValues addObject:[NSNumber numberWithInteger:i]];
        }
    }

    return snpFilterValues;
}


@end
