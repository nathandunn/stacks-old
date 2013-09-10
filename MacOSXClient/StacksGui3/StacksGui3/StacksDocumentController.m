//
//  StacksDocumentController.m
//  StacksGui3
//
//  Created by Nathan Dunn on 8/30/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "StacksDocumentController.h"

@implementation StacksDocumentController

//- (BOOL)validateUserInterfaceItem:(id <NSValidatedUserInterfaceItem>)anItem {
//    NSLog(@"item value %@",anItem);
//    return [super validateUserInterfaceItem:anItem];
//}

- (id)openDocumentWithContentsOfURL:(NSURL *)url display:(BOOL)displayDocument error:(NSError **)outError {
    NSLog(@"using subclass!!");
    [super openDocumentWithContentsOfURL:url display:displayDocument error:outError];
}

- (BOOL)validateUserInterfaceItem:(id <NSValidatedUserInterfaceItem>)anItem {
    NSLog(@"validating UI item in StacksDoc controller %@",anItem) ;
    return [super validateUserInterfaceItem:anItem];
}

- (BOOL)validateMenuItem:(NSMenuItem *)item {
    NSLog(@"validating in in StacksDoc controller menu item %@",item) ;
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
