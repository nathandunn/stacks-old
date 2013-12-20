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

+ (HaplotypeRepository *)sharedInstance
{
    static HaplotypeRepository *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[HaplotypeRepository alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}

- (HaplotypeMO *)insertHaplotype:(NSManagedObjectContext *)context haplotype:(NSString *)haplotype andOrder:(int)order {
    HaplotypeMO *haplotypeMO = [NSEntityDescription insertNewObjectForEntityForName:@"Haplotype" inManagedObjectContext:context];
    haplotypeMO.haplotype = haplotype ;
    haplotypeMO.order = [NSNumber numberWithInt:order];
    return haplotypeMO ;
}
@end