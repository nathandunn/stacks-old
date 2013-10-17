//
//  LocusMO.m
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import "LocusMO.h"
#import "DatumMO.h"
#import "LocusAlleleMO.h"
#import "LocusSnpMO.h"
#import "SampleMO.h"
#import "DepthMO.h"
#import "HaplotypeMO.h"


@implementation LocusMO

@dynamic basePairs;
@dynamic parentCount;
@dynamic chromosome;
@dynamic strand;
@dynamic consensus;
@dynamic length;
@dynamic locusId;
@dynamic marker;
@dynamic ratio;
@dynamic alleles;
@dynamic datums;
@dynamic snps;

@synthesize haplotypeOrder;

- (id)init {
    self = [super init];
    if (self) {
    }

    return self;
}


- (NSInteger) countParents{
    NSInteger count =0 ;
    for(DatumMO *datumMO in self.datums){
        NSString* sampleName = datumMO.sample.name ;
        if([sampleName rangeOfString:@"male"].location!=NSNotFound){
            ++count ;
        }
    }
    return count  ;
}

- (NSInteger) countProgeny{
    NSInteger count =0 ;
//    NSLog(@"number of datums!! %ld",self.datums.count);
    for(DatumMO *datumMO in self.datums){
        NSString* sampleName = datumMO.sample.name ;
//        NSLog(@"sample name %@",sampleName);
        if([sampleName rangeOfString:@"male"].location==NSNotFound){
//            NSLog(@"not found!!") ;
            ++count ;
        }
    }
    return count  ;
}

- (NSAttributedString *)renderChromosome{
    NSString *inputString = @"";
    NSLog(@"self.chromosome %@",self.chromosome);
    if(self.chromosome!=nil && self.chromosome.length>0){
        inputString = [NSString stringWithFormat:@"%@ %@ Mb %@",self.chromosome,self.basePairs,self.strand];
    }
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:inputString];
    return string ;
}

- (NSString *)renderDescription{
    NSUInteger parentCount = [self countParents];
    NSUInteger progenyCount = [self countProgeny];
    NSUInteger snpCount = self.snps.count;
    NSString *chromosomeString = @"";
    if(self.chromosome!=nil && self.chromosome.length>0){
        chromosomeString = [NSString stringWithFormat:@"%@ %@ Mb %@",self.chromosome,self.basePairs,self.strand];
    }
    

    NSString *inputString = [NSString stringWithFormat:@"Parents %ld Prog %ld Snps %ld %@",parentCount,progenyCount,snpCount,chromosomeString];

    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:@"abc123"];
    [string beginEditing];
    [string endEditing];
    return inputString ;
}

- (NSAttributedString *)renderConsensus{
    NSMutableAttributedString *string = [[NSMutableAttributedString alloc] initWithString:self.consensus];
    [string beginEditing];
//        NSNumber *snpIndex;
    NSDictionary *attributes = [NSDictionary dictionaryWithObjectsAndKeys:
            [NSColor blueColor], NSForegroundColorAttributeName,
            [NSColor grayColor], NSBackgroundColorAttributeName,
            [NSFont fontWithName:@"Courier" size:14.0], NSFontAttributeName,
            nil];
    for (LocusSnpMO *snp in self.snps) {
        NSRange selectedRange = NSMakeRange([snp.column unsignedIntegerValue], 1);
        [string setAttributes:attributes range:selectedRange];
    }
    [string endEditing];
    return string ;
}

- (NSUInteger) lookupHaplotypeOrder:(NSString *)haplotype {
    if(haplotypeOrder==nil){
        haplotypeOrder = [[NSMutableDictionary alloc] init];
        [haplotypeOrder setValue:[NSNumber numberWithInt:0] forKey:haplotype];
        return 0 ;
    }

    NSNumber *returnType = [haplotypeOrder objectForKey:haplotype];
    if(returnType==nil){
        NSUInteger maxValue = 0 ;
        for(NSNumber *aValue in haplotypeOrder.allValues){
           if([aValue unsignedIntegerValue]>maxValue) {
               maxValue = [aValue unsignedIntegerValue];
           }
        }
        ++maxValue;
        [haplotypeOrder setValue:[NSNumber numberWithInt:maxValue] forKey:haplotype];
        // actually we get he max value here first .

        return maxValue ;
    }
    else{
        return [returnType unsignedIntegerValue];
    }

}

- (NSUInteger)lookupDepthOrder:(DepthMO *)depthMO{
    // assume that it is there . . .
    NSNumber *order = depthMO.order;
    for(HaplotypeMO *haplotypeMO in depthMO.datum.haplotypes.allObjects){
        if([haplotypeMO.order isEqualToNumber:order]){
            return [self lookupHaplotypeOrder:haplotypeMO.haplotype];
        }
    }

    return 0;
}
@end
