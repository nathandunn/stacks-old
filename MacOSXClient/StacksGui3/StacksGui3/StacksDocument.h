//
//  StacksDocument.h
//  StacksGui3
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface StacksDocument : NSPersistentDocument

@property (strong) NSMutableDictionary *locusViews;
@property(strong) NSString *path;
@property(atomic, retain) NSMutableDictionary *populationLookup;
@property(atomic, retain) NSArray *orderedLocus;
@property(nonatomic, strong) NSSet *loci;

- (id)initWithLocusView:(NSMutableDictionary*)locusViews;
- (id)initWithLoci:(NSSet*)loci ;
- (NSMutableArray *)findPopulations;


- (NSManagedObjectContext *)getContext;

- (NSManagedObjectContext *)getContextForPath:(NSString *)path;

- (NSManagedObjectContext *)getContextForPath:(NSString *)string andName:(NSString *)name;
@end
