//
// Created by ndunn on 2/27/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class GenotypeEntry;

/**
* This is the popup view from each catalog:
* http://genome.uoregon.edu/stacks/catalog.php?db=tut_radtags&id=1
*/
@interface GenotypeView : NSObject

// list of GenotypeEntry's as progeny
@property (strong) NSMutableArray *progeny;

@property (strong) GenotypeEntry *male;
@property (strong) GenotypeEntry *female;

@end