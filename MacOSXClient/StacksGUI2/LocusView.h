//
//  LocusView.h
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>

@class GenotypeView;
@class StacksView;

/**
* This is the table view:
* http://genome.uoregon.edu/stacks/catalog.php?db=tut_radtags&id=1
*/
@interface LocusView : NSObject

@property (strong) NSString *locusId;
//@property (assign) NSString *letters;
@property (strong) NSMutableArray *snps;
@property (strong) NSString *consensus;
@property (strong) NSString *matchingParents;
@property (strong) NSString *progeny;
@property (strong) NSString *marker ;
@property (strong) NSString *ratio ;
@property (strong) NSString *genotypes ;

//@property (strong) GenotypeView *genotypeView;
//@property (strong) NSMutableArray *snpsViews; // SNP in the consensus

// TODO: look in populations.cc
// use PopMap to get the GenotypesView

- (id)initWithId:(NSString *)locusId ;
- (id)initWithId:(NSString *)locusId consensus:(NSString*)consensus;

@end
