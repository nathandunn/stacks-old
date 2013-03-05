//
// Created by NathanDunn on 3/4/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class GenotypeEntry;


@interface DataStubber : NSObject
- (NSMutableArray *)generateSnps;
- (GenotypeEntry *)generateGenotype;

- (NSMutableArray *)generateProgeny:(NSInteger) totalGenotypes;

- (NSString *)generateMarker;
@end