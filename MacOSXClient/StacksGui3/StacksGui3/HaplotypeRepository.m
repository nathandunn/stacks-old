//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "HaplotypeRepository.h"
#import "HaplotypeMO.h"


@implementation HaplotypeRepository {

}
- (HaplotypeMO *)insertHaplotype:(NSManagedObjectContext *)context haplotype:(NSString *)haplotype {
    HaplotypeMO *haplotypeMO = [NSEntityDescription insertNewObjectForEntityForName:@"Haplotype" inManagedObjectContext:context];
    haplotypeMO.haplotype = haplotype ;
    return haplotypeMO ;
}
@end