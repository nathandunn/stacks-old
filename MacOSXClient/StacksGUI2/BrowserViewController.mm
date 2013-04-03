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
#import "LocusView.h"
#import "GenotypeView.h"
#import "GenotypeEntry.h"

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
//    NSLog(@"number of children for class name %@", [item className]);
    NSString *className = [item className];
    if([className isEqualToString:@"StacksDocument"]){
        return ((StacksDocument *) item).childCount;
    }
    else
    if([className isEqualToString:@"LocusView"]){
        LocusView* locusView = (LocusView*) item;
        NSUInteger genotypeCount = locusView.childCount;
//        NSLog(@"count return %ld",genotypeCount);
        
        
        return genotypeCount;
    }
    else{
//        NSLog(@"returing 0 cause I don't know the class %@",className);
        return 0 ;
    }
}

- (id)browser:(NSBrowser *)browser child:(NSInteger)index ofItem:(id)item {
//    NSLog(@"getting child %ld %@",index,item);
//    FileSystemNode *node = (FileSystemNode *)item;
//    return [node.children objectAtIndex:index];
//    NSLog(@"class name %@", [item className]);
    NSString *className = [item className];
    if([className isEqualToString:@"StacksDocument"]){
        return [((StacksDocument *) item) childAtIndex:index];
    }
    else
    if([className isEqualToString:@"LocusView"]){
        // get genoyptes
        return [((LocusView*) item) childAtIndex:index];

    }
    // Not sure what would happen here
//    else
//    if([className isEqualToString:@"GenotypeEntry"]){
//        // get genoyptes
//        GenotypeEntry *genotypeEntry = (GenotypeEntry*) item;
//        StacksView *stacksView = [_stacksLoader loadStacksView:genotypeEntry.name atPath:@"/tmp/stacks_tut/" forTag:genotypeEntry.tagId];
//        return stacksView;
//    }
    else{
        return nil ;
    }
}

- (BOOL)browser:(NSBrowser *)browser isLeafItem:(id)item {
//    NSLog(@"leaf class name %@", [item className]);
    NSString *className = [item className];
    if([className isEqualToString:@"StacksDocument"]){
//        return [((StacksDocument *) item) isLeaf];
        return FALSE ;
    }
    else
    if([className isEqualToString:@"LocusView"]){
        return FALSE ;
    }
    else
    if([className isEqualToString:@"GenotypeEntry"]){
        return TRUE;
    }
    else{
        return TRUE ;
    }
}

- (id)browser:(NSBrowser *)browser objectValueForItem:(id)item {
//    NSLog(@"object for value%@",item);
//    NSLog(@"objectValue for item %@", [item className]);
    NSString *className = [item className];
    if([className isEqualToString:@"StacksDocument"]){
        StacksDocument *doc = ((StacksDocument *) item);
        return [NSString stringWithFormat:@"stacks doc - %@",doc.name] ;
    }
    else
    if([className isEqualToString:@"LocusView"]){
        // get genoyptes
        LocusView *locusView = (LocusView*) item;
        return [NSString stringWithFormat:@"locus - %@",locusView.locusId] ;
    }
    else
    if([className isEqualToString:@"GenotypeEntry"]){
        GenotypeEntry *genotypeEntry = (GenotypeEntry*) item;
//        return genotypeEntry.name ;
        return [NSString stringWithFormat:@"genotype - %@",genotypeEntry.name] ;
        //        return TRUE;
    }
    else{
        return [NSString stringWithFormat:@"other class %@",[item className]];
    }
//    FileSystemNode *node = (FileSystemNode *)item;
//    return node.displayName;
}


- (void)awakeFromNib {

//    [self.browser setDelegate:self];
}


@end
