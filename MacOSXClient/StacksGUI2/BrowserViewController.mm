//
//  BrowserViewController.m
//  StacksGUI2
//
//  Created by Nathan Dunn on 4/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "BrowserViewController.h"
#import "StacksDocument.h"
#import "StacksLoader.h"

@interface BrowserViewController()


@end


@implementation BrowserViewController

@synthesize stacksDocument;



// This method is optional, but makes the code much easier to understand
- (id)rootItemForBrowser:(NSBrowser *)browser {
    // not sure what this would do per se

//    if (_rootNode == nil) {
//        _rootNode = [[FileSystemNode alloc] initWithURL:[NSURL fileURLWithPath:@"/"]];
//    }
//    return _rootNode;
    return self.stacksDocument;
}

- (NSInteger)browser:(NSBrowser *)browser numberOfChildrenOfItem:(id)item {
    NSLog(@"class name %@", [item className]);
    NSString *className = [item className];
    if([className isEqualToString:@"StacksDocument"]){
        return ((StacksDocument *) item).childCount;
    }
    else{
        return 0 ;
    }
}

- (id)browser:(NSBrowser *)browser child:(NSInteger)index ofItem:(id)item {
    NSLog(@"getting child %ld %@",index,item);
//    FileSystemNode *node = (FileSystemNode *)item;
//    return [node.children objectAtIndex:index];
    NSLog(@"class name %@", [item className]);
    NSString *className = [item className];
    if([className isEqualToString:@"StacksDocument"]){
        return [((StacksDocument *) item) childAtIndex:index];
    }
    else{
        return nil ;
    }
}

- (BOOL)browser:(NSBrowser *)browser isLeafItem:(id)item {
    NSLog(@"leaf class name %@", [item className]);
    NSString *className = [item className];
    if([className isEqualToString:@"StacksDocument"]){
        return [((StacksDocument *) item) isLeaf];
    }
    else{
        return TRUE ;
    }
}

- (id)browser:(NSBrowser *)browser objectValueForItem:(id)item {
    NSLog(@"object for value%@",item);
    NSLog(@"objectValue for item %@", [item className]);
    NSString *className = [item className];
    if([className isEqualToString:@"StacksDocument"]){
        StacksDocument *doc = ((StacksDocument *) item);
        return @"a stacks doc";
    }
    else{
        return @"crapola";
    }
//    FileSystemNode *node = (FileSystemNode *)item;
//    return node.displayName;
}


- (void)awakeFromNib {

//    [self.browser setDelegate:self];
}


@end
