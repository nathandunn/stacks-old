//
//  StacksDocument.h
//  StacksGui3
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@class DatumMO;
@class LocusMO;
@class PopulationMO;

@interface StacksDocument : NSPersistentDocument

@property(strong) NSString *datumPath;
@property(strong) NSString *path;
@property(strong) NSString *importPath;
@property(atomic, retain) NSMutableDictionary *populationLookup;

@property(nonatomic, strong) NSSet *loci;
@property(nonatomic, strong) NSSet *populations;


// handle selections
@property(nonatomic, strong) LocusMO *selectedLocus;
@property(nonatomic, strong) PopulationMO *selectedPopulation;
@property(nonatomic, strong) NSArray *selectedDatums;
@property(nonatomic, strong) DatumMO *selectedDatum;

// repositories


@property(nonatomic, retain) NSString *previousStacksName;
@property(nonatomic, retain) NSMutableArray *snpFilterValues;
@property(nonatomic, retain) NSMutableArray *sampleFilterValues;


@property(nonatomic, copy) NSString *name;

// this will be GeneticMap or Population
@property(nonatomic, strong) NSString *type;

@property(nonatomic, strong) NSSet *lociLocations;
@property NSUInteger maxLocation;
@property(strong) NSString *oldPopulationTitle;
@property(nonatomic, strong) NSNumberFormatter *numberFormatter;

- (BOOL)noLociLocations;

- (IBAction)updateSelections:(id)sender;

- (NSArray *) getSnpFilterValues;

@end
