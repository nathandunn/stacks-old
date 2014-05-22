//
//  StacksDocumentController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 8/30/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "StacksDocumentController.h"

@implementation StacksDocumentController

- (id)openDocumentWithContentsOfURL:(NSURL *)url display:(BOOL)displayDocument error:(NSError **)outError {
    return [super openDocumentWithContentsOfURL:url display:displayDocument error:outError];
}

- (BOOL)validateUserInterfaceItem:(id <NSValidatedUserInterfaceItem>)anItem {
    return [super validateUserInterfaceItem:anItem];
}

- (BOOL)validateMenuItem:(NSMenuItem *)item {
//    return [super validateUserInterfaceItem:item];
    if(item.tag==77){
        return YES;
    }
    else{
        return [super validateMenuItem:item];
    }
}


@end
