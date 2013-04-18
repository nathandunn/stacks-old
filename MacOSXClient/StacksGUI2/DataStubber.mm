//
// Created by NathanDunn on 3/4/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "DataStubber.h"
#import "stacks.h"
#import "GenotypeEntry.h"


@implementation DataStubber {

}

- (NSMutableArray *)generateSnps {
    // data Stub!!!
    NSMutableArray *snps = [[NSMutableArray alloc] init];
    while (arc4random_uniform(2) > 0) {
        SNP *snp = new SNP();
        snp->col = 19;
        snp->type = snp_type::snp_type_het;
        snp->rank_1 = 'T';
        snp->rank_2 = 'C';
        snp->lratio = 47.8;
        // may not want to do it this way . . . but we'll see
        [snps addObject:[NSValue valueWithPointer:snp]];
    }
    return snps;
}

/**
* a 95% change of creating
*/
- (GenotypeEntry *)generateGenotype {
    if (arc4random_uniform(100) > 5) {
        GenotypeEntry *genotypeEntry = [[GenotypeEntry alloc] init];
        genotypeEntry.superScript = [self generateLetter];
        genotypeEntry.subScript = [self generateLetter];
        // one has to be valid .  ..
        if(genotypeEntry.subScript==nil && genotypeEntry.superScript==nil){
            genotypeEntry.subScript=@"T";
        }

//        genotypeEntry.sampleId=
//        genotypeEntry.sampleId= 12;
        return genotypeEntry;
    }
    return nil;
}

/**
* 50% change of creating a letter
*/
- (NSString *)generateLetter {
    int number = arc4random_uniform(6);
    switch (number) {
        case 0:
            return @"T";
        case 1:
            return @"C";
        case 2:
            return @"A";
        case 3:
            return @"G";
        default:
            return nil;
    }
}

- (NSMutableArray *)generateProgeny:(NSInteger)totalGenotypes {
    NSMutableArray *progeny = [[NSMutableArray alloc] init];
    for (int i = 0; i < totalGenotypes; i++) {
        GenotypeEntry *genotypeEntry = [self generateGenotype];
        if (genotypeEntry != nil) {
            genotypeEntry.sampleId =i;
            [progeny addObject:genotypeEntry];
        }
    }
    return progeny;
}

- (NSMutableDictionary *)generateGenotypes:(NSInteger)totalGenotypes {
    NSMutableDictionary *progeny = [[NSMutableDictionary alloc] initWithCapacity:totalGenotypes];
    for (int i = 0; i < totalGenotypes; i++) {
        GenotypeEntry *genotypeEntry = [self generateGenotype];
        if (genotypeEntry != nil) {
            genotypeEntry.sampleId =i;
            [progeny setObject:genotypeEntry forKey:[NSString stringWithFormat:@"%ld", genotypeEntry.sampleId] ];
        }
    }
    return progeny;
}

- (NSString *)generateMarker {
    int number = arc4random_uniform(10);
    switch (number) {
        case 0: return @"aa/bb";
        case 1: return @"ab/--";
        case 2: return @"--/ab";
        case 3: return @"aa/ab";
        case 4: return @"ab/aa";
        case 5: return @"ab/ab";
        case 6: return @"ab/ac";
        case 7: return @"ab/cd";
        case 8: return @"ab/cc";
        case 9: return @"cc/ab";
        default:
            return nil;
    }
}
@end