//
//  LocusMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class AlleleMO, GenotypeMO, SnpMO;

@interface LocusMO : NSManagedObject

@property (nonatomic, retain) NSString * consensus;
@property (nonatomic, retain) NSNumber * locusId;
@property (nonatomic, retain) NSString * marker;
@property (nonatomic, retain) NSNumber * length;
@property (nonatomic, retain) NSString * ratio;
@property (nonatomic, retain) NSSet *genotypes;
@property (nonatomic, retain) SnpMO *snps;
@property (nonatomic, retain) AlleleMO *alleles;
@end

@interface LocusMO (CoreDataGeneratedAccessors)

- (void)addGenotypesObject:(GenotypeMO *)value;
- (void)removeGenotypesObject:(GenotypeMO *)value;
- (void)addGenotypes:(NSSet *)values;
- (void)removeGenotypes:(NSSet *)values;

@end
