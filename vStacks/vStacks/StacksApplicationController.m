//
// Created by Nathan Dunn on 9/9/13.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//
//
//


#import "StacksApplicationController.h"
#import "StacksConverter.h"
#import "StacksDocumentController.h"
#import "StacksDocument.h"
#import "ProgressController.h"
#import "StacksAppDelegate.h"
#import "SampleRepository.h"
#import "SampleMO.h"
#import "PopulationRepository.h"
#import "PopulationMO.h"
#import "PopulationService.h"


@implementation StacksApplicationController {


}


@synthesize numberFormatter;

- (id)init {
    self = [super init];
    if (self) {

        numberFormatter = [[NSNumberFormatter alloc] init];
        numberFormatter.numberStyle = NSNumberFormatterNoStyle;
    }

    return self;
}

- (IBAction)importDocument:(id)sender {

    NSLog(@"Importing Document");


    NSOpenPanel *importPanel = [NSOpenPanel openPanel];
    [importPanel setAllowsMultipleSelection:NO];

    [importPanel setCanChooseDirectories:YES];
    [importPanel setCanChooseFiles:NO];

    [importPanel setFloatingPanel:YES];

    NSSize minSize;
    minSize.height = 600;
    minSize.width = 500;


    [importPanel setMinSize:minSize];
    NSInteger result = [importPanel runModal];
    NSLog(@"import result %ld", result);


    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
//    NSInteger result = [panel runModalForDirectory:NSHomeDirectory() file:nil types:nil];
    if (result == NSOKButton) {
        NSString *importPath = [importPanel.directoryURL.path stringByAppendingString:@"/"];
        NSLog(@"import path %@", importPath);
        NSString *importPathName = importPanel.directoryURL.lastPathComponent;


        NSLog(@"import path name %@", importPathName);

        NSFileManager *fileManager = [NSFileManager defaultManager];
        NSArray *files = [fileManager contentsOfDirectoryAtPath:importPath error:nil];
        NSLog(@"# of files %ld", files.count);
        for (id file in files) {
            NSLog(@"file %@", file);
        }

        NSString *fileName = [importPathName stringByAppendingString:@".stacks"];
        NSString *extension = [fileName pathExtension];
        if ([extension isNotEqualTo:@"stacks"]) {
            fileName = [fileName stringByAppendingString:@".stacks"];
            NSLog(@"has correct filename %i", [[fileName exposedBindings] isEqualTo:@"stacks"]);
            NSLog(@"filename: %@", fileName);
        }
        NSString *savedStacksDocumentPath = [importPath stringByAppendingFormat:@"%@", fileName];
        NSLog(@"stacks doc path %@", savedStacksDocumentPath);
        NSString *resultFileName ;

        if ([[NSFileManager defaultManager] fileExistsAtPath:savedStacksDocumentPath isDirectory:NULL]) {
            NSAlert *alertReplace = [[NSAlert alloc] init];
            [alertReplace setMessageText:[NSString stringWithFormat:@"Replace stacks file or create new%@?", savedStacksDocumentPath]];
            [alertReplace addButtonWithTitle:@"Cancel"];
            [alertReplace addButtonWithTitle:@"Replace"];
            [alertReplace addButtonWithTitle:@"Create New"];

            NSInteger runAlertReplace = [alertReplace runModal];
            if (runAlertReplace == NSAlertFirstButtonReturn) {
                NSLog(@"Cancelling");
                return;
            }
            if (runAlertReplace == NSAlertSecondButtonReturn) {
                BOOL fileRemoved = [[NSFileManager defaultManager] removeItemAtPath:savedStacksDocumentPath error:NULL];
                NSLog(@"file removed %i", fileRemoved);
            }
        }

        int i = 1;
        while ([[NSFileManager defaultManager] fileExistsAtPath:savedStacksDocumentPath isDirectory:NULL]) {
            resultFileName = [importPathName stringByAppendingFormat:@"%i.stacks", i];
            savedStacksDocumentPath = [importPath stringByAppendingFormat:@"%@", resultFileName];
            ++i;
        }
//        fileName = resultFileName ;


        NSAlert *alertCancel = [[NSAlert alloc] init];
        [alertCancel setMessageText:[NSString stringWithFormat:@"Create stacks file %@?", savedStacksDocumentPath]];
        [alertCancel addButtonWithTitle:@"Cancel"];
        [alertCancel addButtonWithTitle:@"OK"];

        NSInteger runAlert = [alertCancel runModal];
        if (runAlert == NSAlertFirstButtonReturn) {
            NSLog(@"Cancelling");
            return;
        }

//        if (fileExistsAtPath) {
//            BOOL fileRemoved = [[NSFileManager defaultManager] removeItemAtPath:savedStacksDocumentPath error:NULL];
//            NSLog(@"file removed %i", fileRemoved);
//        }


        ProgressController *progressController = [[ProgressController alloc] init];

        progressController.stacksConverter = stacksConverter;
        [progressController showWindow:[NSApp mainWindow]];

//        dispatch_async(dispatch_get_main_queue(),^ {
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{

            NSLog(@"loadding progress!!! in thread");
//            StacksDocument *newDocument = [stacksConverter loadLociAndGenotypes:[panel.directoryURL.path stringByAppendingString:@"/"] progressWindow:progressController];
            StacksDocument *newDocument = [stacksConverter loadLociAndGenotypes:savedStacksDocumentPath progressWindow:progressController importPath:importPath];
            newDocument.path = savedStacksDocumentPath;
            [newDocument.managedObjectContext save:nil];
            if (newDocument != nil) {
                [progressController close];
                NSLog(@"LOADED progress!!! in thread");
//        [NSApp stopModal];
                NSLog(@"trying to open");

                [[StacksDocumentController sharedDocumentController] openDocumentWithContentsOfURL:[NSURL fileURLWithPath:savedStacksDocumentPath] display:YES completionHandler:^(NSDocument *doc, BOOL documentWasAlreadyOpened, NSError *error) {
                    if (error != nil) {
                        NSLog(@"error3 %@", error);
                    }
                }];


                for (StacksDocument *stacksDocument in [[StacksDocumentController sharedDocumentController] documents]) {
                    NSLog(@"stacks doc: %@", stacksDocument.path);
                    if (stacksDocument.path == NULL) {
                        [stacksDocument close];
//                [[StacksDocumentController sharedDocumentController] perform]
                    }
                }

            }
            else {
                [progressController close];
                NSLog(@"must have been cancelled");

            }
        }
        );
//        });
//        [stacksConverter loadLociAndGenotypes:[panel.directoryURL.path stringByAppendingString:@"/"] progressWindow:progressController];



    }
    else {
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
//- (void)startProgressPanel:(NSString *)message {
//    [loadProgress displayIfNeeded];
//    [loadProgress setIndeterminate:false];
//    [loadProgress setDisplayedWhenStopped:false];
//    [loadProgress setNeedsDisplay:true];

//    [self beginSheetModalForWindow:[ progressPanel window]  modalDelegate:self  didEndSelector:@ selector( alertEnded:code:context:)  contextInfo:NULL];


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
//}

//- (void)progressDidEnd:(NSWindow *)panel returnCode:(int)returnCode contextInfo:(void *)context {
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
//}


// TODO: change for file, etc. etc. etc.
- (IBAction)applyPopmap:(id)sender {

    StacksDocument *stacksDocument = [[StacksDocumentController sharedDocumentController] currentDocument];
    NSLog(@"stacksDocu?  %@", stacksDocument);
    if (stacksDocument == nil) {
        NSAlert *alert = [NSAlert alertWithMessageText:@"Error applying Population Map"
                                         defaultButton:@"OK"
                                       alternateButton:nil
                                           otherButton:nil
                             informativeTextWithFormat:@"Must have an open Stacks document in order to apply a Population Map"];
        [alert runModal];
        return;
    }

    NSOpenPanel *importPanel = [NSOpenPanel openPanel];
    [importPanel setAllowsMultipleSelection:NO];
    [importPanel setCanChooseDirectories:NO];
    [importPanel setCanChooseFiles:YES];
    [importPanel setFloatingPanel:YES];
    NSSize minSize;
    minSize.height = 600;
    minSize.width = 500;


    [importPanel setMinSize:minSize];
    NSInteger result = [importPanel runModal];


    if (result == NSOKButton) {
        NSString *fileUrlString = [NSString stringWithFormat:@"%@", [importPanel.URL path]];
        NSLog(@"import file %@", fileUrlString);

        NSString *errorCondition = [[PopulationService sharedInstance] validatePopmap:importPanel.URL];
        if (errorCondition != nil) {
            NSLog(@"Bad Popmap %@", errorCondition);

            NSAlert *alert = [NSAlert alertWithMessageText:[NSString stringWithFormat:@"Unable to apply Population Map %@", importPanel.URL.path]
                                             defaultButton:@"OK"
                                           alternateButton:nil
                                               otherButton:nil
                                 informativeTextWithFormat:errorCondition];

            [alert runModal];

            return;
        }



        // remove the old populations from the samples
        NSArray *sampleArray = [[SampleRepository sharedInstance] getAllSamples:stacksDocument.managedObjectContext];
        for (SampleMO *sampleMO in sampleArray) {
            sampleMO.population = nil ;
        }

        // remove the old populations
        NSArray *allPopulation = [[PopulationRepository sharedInstance] getAllPopulations:stacksDocument.managedObjectContext];
        for (PopulationMO *populationMO in allPopulation) {
            [stacksDocument.managedObjectContext deleteObject:populationMO];
        }

        [stacksDocument.populationLookup removeAllObjects];
        stacksDocument.populationLookup = nil ;

        // TODO: validate deletions!!!


        // TODO: refactor to take a filename
        StacksConverter *stacksConverter = [[StacksConverter alloc] init];
        [stacksConverter addPopulationsToDocument:stacksDocument forPath:fileUrlString];
        NSLog(@"population lookup size: %ld", stacksDocument.populationLookup.count);

        // TODO: reapply to the samples
        for (SampleMO *sampleMO in sampleArray) {
            NSString *populationId = [stacksDocument.populationLookup objectForKey:sampleMO.name];

            if (populationId != nil) {
                // lets get the population . . can use lookup, but this is usually pretty small
                NSLog(@"tyring to populate popid %@", populationId);
                for (PopulationMO *populationMO in stacksDocument.populations) {
//                    NSNumber *endNumber = [numberFormatter numberFromString:populationId];
//                    NSLog(@"comparing to %@ vs %@", populationMO.populationId, endNumber);
                    if ([populationMO.name isEqualToString:populationId]) {
                        NSLog(@"FOuND population ID %@", populationMO.populationId);
                        [populationMO addSamplesObject:sampleMO];
                    }
                }
            }
        }


        NSError *saveError = nil ;
        [stacksDocument.managedObjectContext save:&saveError];
        if (saveError != nil) {
            NSLog(@"save error %@", saveError);
        }


        NSURL *url = stacksDocument.fileURL;
        NSError *openError = nil ;
        [[self currentDocument] close];
        [self reopenDocumentForURL:url withContentsOfURL:url error:&openError];
//        [self openDocumentWithContentsOfURL: display:<#(BOOL)displayDocument#> error:<#(NSError **)outError#>:stacksDocument];

        // TODO: reapply to the  UI


    }

}


- (void)stopProgressPanel {

//    [progressPanel orderOut: self];
//    [NSApp endSheet: progressPanel returnCode: 0];
}

- (IBAction)cancelAction:(id)sender {
    NSLog(@"in stacksappcontroller . .. cancelling action");
//    [progressPanel orderOut: self];
//    [NSApp endSheet: progressPanel returnCode: 1];
}

//
- (void)someMethodDidEnd:(NSAlert *)alert returnCode:(int)returnCode contextInfo:(void *)contextInfo {
    if (returnCode == NSAlertFirstButtonReturn) {
        // Do something
    }
}

- (void)provideFeedback:(id)sender {
    NSString *recipients = @"mailto:jcatchen@uoregon.edu?cc=ndunn@me.com&subject=vStacks Feedback 0.8.4";
    NSString *body = @"&body=Feedback for vStacks";
    NSString *email = [NSString stringWithFormat:@"%@%@", recipients, body];

    email = [email stringByAddingPercentEscapesUsingEncoding:NSUTF8StringEncoding];

    [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:email]];
}

- (void)helpStacks:(id)sender {
    NSString *url = @"http://creskolab.uoregon.edu/stacks/";
    [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:url]];
}

- (void)helpVStacks:(id)sender {
    NSString *url = @"http://creskolab.uoregon.edu/stacks/manual/#vstacks";
    [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:url]];
}

- (void)license:(id)sender {
    NSString *url = @"http://creskolab.uoregon.edu/stacks/vstacks/apple-license.php";
    [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:url]];
}

@end