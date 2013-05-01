//
// Created by ndunn on 2/27/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class StackEntryMO;

/**
* This is just a raw read formatted:
* http://genome.uoregon.edu/stacks/tag.php?db=tut_radtags&batch_id=1&sample_id=7&tag_id=408
*/
@interface StacksView : NSObject

// a list of StacksEntry's
@property (strong) NSMutableArray *stackEntries;
@property (strong) StackEntryMO *consensus;
@property (strong) StackEntryMO *model;
@property (strong) StackEntryMO *reference;

// index of snps in the consensus
@property (strong) NSMutableArray *snps;


- (NSInteger)rowsNeeded;
@end