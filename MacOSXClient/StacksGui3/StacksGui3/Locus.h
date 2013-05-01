//
//  Locus.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/1/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class Allele, Genotype, Snp;

@interface Locus : NSManagedObject

@property (nonatomic, retain) NSString * consensus;
@property (nonatomic, retain) NSNumber * locusId;
@property (nonatomic, retain) NSString * marker;
@property (nonatomic, retain) NSNumber * length;
@property (nonatomic, retain) NSString * ratio;
@property (nonatomic, retain) NSSet *genotypes;
@property (nonatomic, retain) Snp *snps;
@property (nonatomic, retain) Allele *alleles;
@end

@interface Locus (CoreDataGeneratedAccessors)

- (void)addGenotypesObject:(Genotype *)value;
- (void)removeGenotypesObject:(Genotype *)value;
- (void)addGenotypes:(NSSet *)values;
- (void)removeGenotypes:(NSSet *)values;

@end
