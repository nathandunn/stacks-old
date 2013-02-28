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

@interface LocusView : NSObject

@property (strong) NSString *locusId;
//@property (assign) NSString *letters;
@property (strong) NSString *snp;
@property (strong) NSString *consensus;
@property (strong) NSString *matchingParents;
@property (strong) NSString *progeny;
@property (strong) NSString *marker ;
@property (strong) NSString *ratio ;
@property (strong) NSString *genotypes ;

@property (strong) GenotypeView *genotypeView;
@property (strong) StacksView *stacksView;

- (id)initWithId:(NSString *)locusId consensus:(NSString*)consensus;

@end
