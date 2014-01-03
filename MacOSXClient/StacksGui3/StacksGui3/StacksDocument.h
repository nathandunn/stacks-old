//
//  StacksDocument.h
//  StacksGui3
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@class DatumMO;
@class LocusMO;
@class PopulationMO;
//@class DatumRepository;
//@class PopulationRepository;
//@class LocusRepository;

@interface StacksDocument : NSPersistentDocument

@property(strong) NSString *path;
@property(atomic, retain) NSMutableDictionary *populationLookup;

@property(nonatomic, copy) NSSet *loci;
@property(nonatomic, copy) NSSet *populations;


// handle selections
@property(nonatomic, strong) LocusMO *selectedLocus;
@property(nonatomic, strong) PopulationMO *selectedPopulation;
@property(nonatomic, copy) NSArray *selectedDatums;
@property(nonatomic, strong) DatumMO *selectedDatum;

// repositories
//@property(nonatomic, strong) DatumRepository *datumRepository;
//@property(nonatomic, strong) LocusRepository *locusRepository;
//@property(nonatomic, strong) PopulationRepository *populationRepository;


@property(nonatomic, retain) NSString *previousStacksName;
@property(nonatomic, copy) NSMutableArray *snpFilterValues;

//- (NSMutableArray *)findPopulations;


@property(nonatomic, copy) NSString *name;

// this will be GeneticMap or Population
@property(nonatomic, strong) NSString *type;

@property(nonatomic, copy) NSSet *lociLocations;
@property NSUInteger maxLocation;

- (BOOL)noLociLocations;

- (IBAction)updateSelections:(id)sender;

- (NSArray *) getSnpFilterValues;

@end
