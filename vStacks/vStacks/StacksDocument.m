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
#import "SampleRepository.h"
#import "SampleMO.h"
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
@property(weak) IBOutlet NSPopUpButton *maxSamplesPopupButton;
@property(weak) IBOutlet WebView *stacksWebView;
@property(weak) IBOutlet WebView *datumWebView;


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
@synthesize sampleFilterValues;
@synthesize maxSnpPopupButton;
@synthesize maxSamplesPopupButton;
@synthesize stacksWebView;
@synthesize datumWebView;
@synthesize numberFormatter;

@synthesize importPath;
@synthesize path;
@synthesize oldPopulationTitle;
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

        numberFormatter = [[NSNumberFormatter alloc] init];
        numberFormatter.numberStyle = NSNumberFormatterNoStyle;
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
    [maxSamplesPopupButton selectItemAtIndex:[self getSampleFilterValues].count - 1];


//    [filteredLoci setFormatter:numberFormatter];
    [[filteredLoci cell] setFormatter:numberFormatter];

    [self updateStacksView];
}

- (void)updateStacksView {

    if (self.selectedDatum != nil) {

        StackEntryDatumMO *stackEntryDatumMO = [[DatumRepository sharedInstance] getStackEntryDatum:self.managedObjectContext datum:self.selectedDatum];
        if (stackEntryDatumMO != nil && stackEntryDatumMO.stackData != nil) {

            NSData *jsonData = [stackEntryDatumMO.stackData gunzippedData];
            StacksEntryDatumRenderer *stacksEntryDatumRenderer = [[StacksEntryDatumRenderer alloc] init];
            NSString *html = [stacksEntryDatumRenderer renderHtmlForData:jsonData datumSnps:self.selectedDatum.snpData locusSnps:self.selectedLocus.snpData];

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
        NSLog(@"hit DONE, setting button to edit");
        editPopulationButton.title = @"Edit";

        for (id item in populationSelector.itemArray) {
            NSLog(@"item in there: %@", item);
        }


//        PopulationMO *populationMO = [[PopulationRepository sharedInstance] getPopulation:self.managedObjectContext name:oldPopulationTitle];
        PopulationMO *populationMO = [[PopulationRepository sharedInstance] getPopulation:self.managedObjectContext name:oldPopulationTitle];
        NSLog(@"popMO: %@", populationMO);
        if (populationMO != nil && ![populationNameField.stringValue isEqualToString:oldPopulationTitle]) {
            populationMO.name = populationNameField.stringValue;
            NSError *error;
            [self.managedObjectContext save:&error];
            if (error) {
                NSLog(@"error");
                return;
            }
        }

        [populationSelector setHidden:false];
        [populationNameField setHidden:true];

        [populationSelector selectItemAtIndex:previousSelectedItem];
        NSLog(@"selected item index %ld", populationSelector.indexOfSelectedItem);
    }
    else {
        if (populationSelector.indexOfSelectedItem == 0) {
            return;
        }
        NSLog(@"setting to DONE");
        previousSelectedItem = populationSelector.indexOfSelectedItem;
        oldPopulationTitle = populationSelector.titleOfSelectedItem;
        editPopulationButton.title = @"Done";
        [populationSelector setHidden:true];
        [populationNameField setHidden:false];

        populationNameField.stringValue = oldPopulationTitle;

    }

    editingPopulation = !editingPopulation;
}


- (void)tableViewSelectionDidChange:(NSNotification *)aNotification {

    self.selectedLocus = [self findSelectedLocus];
    self.selectedPopulation = [self findSelectedPopulation];

    if (self.selectedLocus != nil) {
        NSLog(@"getting selected locus %@", self.selectedLocus.locusId);
//        NSLog(@"getting selected population %@", self.selectedPopulation.name);
        if (self.selectedPopulation != nil) {
            self.selectedDatums = [[DatumRepository sharedInstance] getDatumsOrdered:self.managedObjectContext locus:self.selectedLocus.locusId andPopulation:self.selectedPopulation];
        }
        else {
            self.selectedDatums = [[DatumRepository sharedInstance] getDatumsOrdered:self.managedObjectContext locus:self.selectedLocus.locusId];
        }
        NSLog(@"got selected Datums: %ld", self.selectedDatums.count);
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
    [self updateDatumView];
    [self updateStacksView];

}

- (void)updateDatumView {
    NSString *cssPath = [[NSBundle mainBundle] pathForResource:@"stacks" ofType:@"css"];
    NSString *cssString = [NSString stringWithContentsOfFile:cssPath encoding:NSUTF8StringEncoding error:NULL];
    NSMutableString *returnHTML = [NSMutableString stringWithFormat:@"<style type='text/css'>%@</style>", cssString];

    NSArray *allPopulations = [[PopulationRepository sharedInstance] getAllPopulations:self.managedObjectContext];

    if (self.selectedPopulation) {
        [returnHTML appendFormat:@"<div class='datum-pop'><div class='population-header'>%@</div>", self.selectedPopulation.annotatedName];
        if (self.selectedDatums.count > 0) {
            for (DatumMO *datum in self.selectedDatums.reverseObjectEnumerator) {
                [returnHTML appendString:[self renderDatumHtml:datum]];
            }
        }
        else {
            [returnHTML appendFormat:@"<div class='none'>None</div>"];
        }
        [returnHTML appendFormat:@"</div>"];
    }
    else if (self.selectedPopulation == nil && allPopulations != nil && allPopulations.count > 0) {
        for (PopulationMO *populationMO in allPopulations) {
            NSArray *datums = [[DatumRepository sharedInstance] getDatums:self.managedObjectContext locus:self.selectedLocus.locusId andPopulation:populationMO];
            if (datums.count > 0) {
                [returnHTML appendFormat:@"<div class='datum-pop'><div class='population-header'>%@</div>", populationMO.annotatedName];
                for (DatumMO *datum in datums.reverseObjectEnumerator) {
                    [returnHTML appendString:[self renderDatumHtml:datum]];
                }
                [returnHTML appendFormat:@"</div>"];
            }
        }
    }
            // if no population!
    else {
        [returnHTML appendFormat:@"<div class='datum-pop'><div class='population-header'>Population - Unspecified</div>"];
        for (DatumMO *datum in self.selectedDatums.reverseObjectEnumerator) {
            [returnHTML appendString:[self renderDatumHtml:datum]];
        }
        [returnHTML appendFormat:@"</div>"];
    }

//    NSLog(@"return HTML: %@",returnHTML);

    [[datumWebView mainFrame] loadHTMLString:returnHTML baseURL:[[NSBundle mainBundle] bundleURL]];
}

- (NSString *)renderDatumHtml:(DatumMO *)datum {
    NSMutableString *returnHTML = [NSMutableString string];
    NSString *nameString = [datum renderNameHtml];
//    NSString *haploytpeString = [datum renderHaplotypeHtml] ;
//    NSString *haploytpeString = [[datum renderHaplotypes] string];
    NSString *haploytpeString = [datum renderHaplotypeHtml];
//    NSString *depthString = [[datum renderDepths] string];
    NSString *depthString = [datum renderDepthHtml];
    NSString *datumIndex = [NSString stringWithFormat:@"%@:%@", datum.tagId, datum.sampleId];
    NSString *selectedClass = ([datum isEqualTo:self.selectedDatum]) ? @" selected-datum" : @"";
    [returnHTML appendFormat:@"<div class='population%ld datum%@'>", self.selectedPopulation.populationId.longValue, selectedClass];
    [returnHTML appendFormat:@"<a href='%@'>%@</a><br/>", datumIndex, nameString];
    [returnHTML appendFormat:@"<a href='%@'>%@</a><br/>", datumIndex, haploytpeString];
    [returnHTML appendFormat:@"<a href='%@'>%@</a><br/>", datumIndex, depthString];
    [returnHTML appendFormat:@"</div>"];
    return returnHTML;
}

- (void)webView:(WebView *)sender decidePolicyForNavigationAction:(NSDictionary *)actionInformation
        request:(NSURLRequest *)request frame:(WebFrame *)frame decisionListener:(id)listener {
    NSLog(@"handling click policy!!! %@", [[request URL] lastPathComponent]);

    [listener use];
    NSArray *pathComponents = [[[request URL] lastPathComponent] componentsSeparatedByString:@":"];
    if (pathComponents.count != 2) {
        NSLog(@"path path for URL %@", request.URL);
        return;
    }
    NSUInteger locusId = [numberFormatter numberFromString:[pathComponents objectAtIndex:0]].unsignedIntegerValue;
    NSUInteger sampleId = [numberFormatter numberFromString:[pathComponents objectAtIndex:1]].unsignedIntegerValue;
    DatumMO *datumMO = [[DatumRepository sharedInstance] getDatum:self.managedObjectContext locusId:locusId andSampleId:sampleId];
    if (![datumMO isEqualTo:self.selectedDatum]) {
        self.selectedDatum = datumMO;
        [self updateStacksView];
        [self updateDatumView];
    }


//    NSString *host = [[request URL] host];
//    if ([host hasSuffix:@"company.com"])
//        [listener ignore];
//    else
//        [listener use];
}

- (PopulationMO *)findSelectedPopulation {
    NSInteger selectedRow = [self.populationSelector indexOfSelectedItem];
    if (selectedRow > 0) {
        return [[PopulationRepository sharedInstance] getPopulation:self.managedObjectContext byIndexSortedByName:selectedRow - 1];
    }
    return nil;
}

- (LocusMO *)findSelectedLocus {
    NSInteger selectedRowIndex = [self.locusTableView selectedRow];
    if (selectedRowIndex < 0) {
        return nil;
    }
    NSTableCellView *selectedRow = [self.locusTableView viewAtColumn:0 row:selectedRowIndex makeIfNecessary:YES];
    NSArray *subviews = selectedRow.subviews;
    NSInteger locusId = -1;
    for (NSTextField *subview in subviews) {
        if ([subview.identifier isEqualToString:@"LocusId"]) {
            locusId = [numberFormatter numberFromString:subview.stringValue].integerValue;
        }
    }
    if (locusId >= 0) {
        // id starts at 1 + row  . . . I hope this is always true
        return [[LocusRepository sharedInstance] getLocus:self.managedObjectContext forId:locusId];
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
    NSString *url = @"http://creskolab.uoregon.edu/stacks";
    [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:url]];
}

- (void)provideFeedback:(id)sender {
    NSString *recipients = @"mailto:jcatchen@uoregon.edu?cc=ndunn@uoregon.edu&subject=vStacks Feedback";
    NSString *body = @"&body=Feedback for vStacks";
    NSString *email = [NSString stringWithFormat:@"%@%@", recipients, body];

    email = [email stringByAddingPercentEscapesUsingEncoding:NSUTF8StringEncoding];

    [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:email]];
}

- (NSArray *)getSnpFilterValues {

    NSUInteger maxLocusSnps = [self getMaxLocusSnps];

    if (snpFilterValues == nil) {
        snpFilterValues = [NSMutableArray array];
        for (int i = 0; i < maxLocusSnps + 1; i++) {
            [snpFilterValues addObject:[NSNumber numberWithInteger:i]];
        }
    }

    return snpFilterValues;
}

- (NSArray *)getSampleFilterValues {

    NSUInteger maxLocusSamples = [self getMaxLocusSamples];

    if (sampleFilterValues == nil) {
        sampleFilterValues = [NSMutableArray array];
        for (int i = 1; i < maxLocusSamples + 1; i++) {
            [sampleFilterValues addObject:[NSNumber numberWithInteger:i]];
        }
    }

    return sampleFilterValues;
}

- (NSUInteger)getMaxLocusSamples {
    NSArray *allLocusArray = [[LocusRepository sharedInstance] getAllLoci:self.managedObjectContext];

    NSUInteger maxLocusSamples = 0;
    for (LocusMO *locusMO in  allLocusArray) {

        if (locusMO.progenyCount.unsignedIntegerValue > maxLocusSamples) {
            maxLocusSamples = locusMO.progenyCount.unsignedIntegerValue;
        }
    }

    return maxLocusSamples;
}

- (NSUInteger)getMaxLocusSnps {

    NSArray *allLocusArray = [[LocusRepository sharedInstance] getAllLoci:self.managedObjectContext];

    NSUInteger maxLocusSnps = 0;
    for (LocusMO *locusMO in  allLocusArray) {
        NSArray *snps = [NSJSONSerialization JSONObjectWithData:locusMO.snpData options:kNilOptions error:nil];
        if (snps.count > maxLocusSnps) {
            maxLocusSnps = snps.count;
        }
    }

    return maxLocusSnps;
}


@end
