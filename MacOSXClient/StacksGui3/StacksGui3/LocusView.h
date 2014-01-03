//
//  LocusView.h
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>

@class GenotypeView;
//@class StacksView;
@class GenotypeEntry;

/**
* This is the table view:
* http://genome.uoregon.edu/stacks/catalog.php?db=tut_radtags&id=1
*/
@interface LocusView : NSObject

@property NSInteger locusId;
@property (copy) NSMutableArray *snps;
@property (copy) NSString *consensus;
@property (copy) NSString *marker ;
@property (copy) NSString *ratio ;


@property (copy) NSMutableDictionary *genotypes;

// use PopMap to get the GenotypesView

@property(nonatomic) int depth;

- (id)initWithId:(NSInteger)locusId ;

- (NSInteger) matchingParents;

//- (BOOL)hasMale;
//- (BOOL)hasFemale;
//- (NSInteger)count;
@end
