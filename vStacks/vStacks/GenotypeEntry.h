//
// Created by ndunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>


@interface GenotypeEntry : NSObject

@property NSString *name ;
@property NSInteger sampleId;

// TODO: remove
@property NSString *superScript;
@property NSString *subScript;

// replaced with this:
@property (retain) NSMutableArray *haplotypes;
@property (retain) NSMutableArray *depths;


// link info
//@property NSUInteger  sampleId;

// will use the selected locus
@property NSInteger  tagId;

- (NSString *)render;

- (id)initWithCapcity:(int)numHaplotypes;

- (NSAttributedString *)renderName;
- (NSAttributedString *)renderHaplotypes;
- (NSAttributedString *)renderDepths;
@end