//
// Created by Nathan Dunn on 9/9/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StacksApplicationController.h"
#import "StacksConverter.h"
#import "StacksDocumentController.h"
#import "StacksDocument.h"
#import "ProgressController.h"
#import "StacksAppDelegate.h"


@implementation StacksApplicationController {


}



//- (id)init {
//    self = [super init];
//    if (self) {
//
//    }
//
//    return self;
//}

//- (BOOL)validateUserInterfaceItem:(id <NSValidatedUserInterfaceItem>)anItem {
////    NSLog(@"validating UI item in App Controller %@",anItem) ;
//    return [super validateUserInterfaceItem:anItem];
//}

//- (BOOL)validateMenuItem:(NSMenuItem *)item {
//    NSLog(@"validating in App Controller menu item %@",item) ;
////    return [super validateUserInterfaceItem:item];
//    if(item.tag==77){
//        NSLog(@"should be returning true!");
//        return YES;
//    }
//    else{
//        return [super validateMenuItem:item];
//    }
//}

- (IBAction)importDocument:(id)sender {
    NSLog(@"Importing doc");
    NSOpenPanel *panel = [NSOpenPanel openPanel];
    [panel setAllowsMultipleSelection:NO];
    [panel setCanChooseDirectories:YES];
    [panel setCanChooseFiles:NO];
    [panel setFloatingPanel:YES];
    NSSize minSize;
    minSize.height = 600;
    minSize.width = 500;

    [panel setMinSize:minSize];
    NSInteger result = [panel runModal];

    StacksConverter *stacksConverter = [[StacksConverter alloc] init];
//    NSInteger result = [panel runModalForDirectory:NSHomeDirectory() file:nil types:nil];
    if (result == NSOKButton) {
        NSLog(@"ok !!");
        NSLog(@"directory URL: %@", panel.directoryURL.path);
        NSString *stacksDocumentPath = [stacksConverter generateFilePathForUrl:panel.directoryURL];
        BOOL fileRemoved = [[NSFileManager defaultManager] removeItemAtPath:stacksDocumentPath error:NULL];
        NSLog(@"file removed %i", fileRemoved);

        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{
            NSLog(@"loadding progress!!! in thread");
            ProgressController *progressController = [[ProgressController alloc] init];
            progressController.stacksConverter = stacksConverter;
            [progressController showWindow:[NSApp mainWindow]];
            if ([stacksConverter loadLociAndGenotypes:[panel.directoryURL.path stringByAppendingString:@"/"] progressWindow:progressController] != nil) {
                [progressController close];
                NSLog(@"LOADED progress!!! in thread");
//        [NSApp stopModal];
                NSLog(@"trying to open");
                StacksDocument *opendDoc = [[StacksDocumentController sharedDocumentController] openDocumentWithContentsOfURL:[NSURL fileURLWithPath:stacksDocumentPath] display:YES error:NULL];
                opendDoc.path = stacksDocumentPath;
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
                NSLog(@"must have been cancelled") ;

            }
        });
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
- (void)startProgressPanel:(NSString *)message {
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
}

- (void)progressDidEnd:(NSWindow *)panel returnCode:(int)returnCode contextInfo:(void *)context {
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


- (void)stopProgressPanel {

//    [progressPanel orderOut: self];
//    [NSApp endSheet: progressPanel returnCode: 0];
}

- (IBAction)cancelAction:(id)sender {

//    [progressPanel orderOut: self];
//    [NSApp endSheet: progressPanel returnCode: 1];
}

//
- (void)someMethodDidEnd:(NSAlert *)alert returnCode:(int)returnCode contextInfo:(void *)contextInfo {
    if (returnCode == NSAlertFirstButtonReturn) {
        // Do something
    }
}

@end