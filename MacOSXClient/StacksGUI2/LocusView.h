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

@property (strong) NSString *locusId;
@property (strong) NSMutableArray *snps;
@property (strong) NSString *consensus;
@property (strong) NSString *marker ;
@property (strong) NSString *ratio ;


@property (strong) NSMutableDictionary *genotypes;

// use PopMap to get the GenotypesView

@property(nonatomic) int depth;

- (id)initWithId:(NSString *)locusId ;

- (NSInteger) matchingParents;

- (BOOL)hasMale;
- (BOOL)hasFemale;

- (NSInteger)count;
@end
