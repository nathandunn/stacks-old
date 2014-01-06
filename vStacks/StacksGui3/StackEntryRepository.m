//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StackEntryRepository.h"
#import "StackEntryMO.h"
#import "ConsensusStackEntryMO.h"
#import "DatumMO.h"
#import "ModelStackEntryMO.h"
#import "ReferenceStackEntryMO.h"
#import "DatumSnpMO.h"
#import "LocusMO.h"


NSString *generateSnpString(NSSet *set);

@implementation StackEntryRepository {

}
- (StackEntryMO *)insertStackEntry:(NSManagedObjectContext *)context entryId:(NSNumber *)entryId relationship:(NSString *)relationship block:(NSString *)block sequenceId:(NSString *)sequenceId sequence:(NSString *)sequence consensus:(NSString*)consensus datum:(DatumMO *)datum {
    StackEntryMO *stackEntryMO = [NSEntityDescription insertNewObjectForEntityForName:@"StackEntry" inManagedObjectContext:context];
    stackEntryMO.entryId = entryId;
    stackEntryMO.relationship = relationship;
    stackEntryMO.block = block ;
    stackEntryMO.sequenceId = sequenceId ;
    stackEntryMO.sequence = sequence ;
    stackEntryMO.datum = datum ;
    stackEntryMO.consensus = consensus;

    // TODO: convert 
//    stackEntryMO.datumSnps = generateSnpString(datum.snps);
//    stackEntryMO.locusSnps = generateSnpString(datum.locus.snps);

    return stackEntryMO ;
}

- (ConsensusStackEntryMO *)insertConsensusStackEntry:(NSManagedObjectContext *)context block:(NSString*)block sequenceId:(NSString*)sequenceId sequence:(NSString*)sequence datum:(DatumMO *)datum {
    ConsensusStackEntryMO *stackEntryMO = [NSEntityDescription insertNewObjectForEntityForName:@"ConsensusStackEntry" inManagedObjectContext:context];
    stackEntryMO.relationship = @"consensus" ;
    stackEntryMO.block = block ;
    stackEntryMO.sequenceId = sequenceId ;
    stackEntryMO.sequence = sequence ;
    stackEntryMO.datum = datum ;

    // TODO: convert
//    stackEntryMO.datumSnps = generateSnpString(datum.snps);
//    stackEntryMO.locusSnps = generateSnpString(datum.locus.snps);

    return stackEntryMO ;
}

- (ModelStackEntryMO *)insertModelStackEntry:(NSManagedObjectContext *)context block:(NSString *)block sequenceId:(NSString *)sequenceId sequence:(id)sequence datum:(DatumMO *)datum {
    ModelStackEntryMO *stackEntryMO = [NSEntityDescription insertNewObjectForEntityForName:@"ModelStackEntry" inManagedObjectContext:context];
    stackEntryMO.relationship = @"model" ;
    stackEntryMO.block = block ;
    stackEntryMO.sequenceId = sequenceId ;
    stackEntryMO.sequence = sequence ;
    stackEntryMO.datum = datum ;

    // TODO: convert 
//    stackEntryMO.datumSnps = generateSnpString(datum.snps);
//    stackEntryMO.locusSnps = generateSnpString(datum.locus.snps);

    return stackEntryMO ;
}



- (ReferenceStackEntryMO *)insertReferenceStackEntry:(NSManagedObjectContext *)context sequence:(NSString *)sequence datum:(DatumMO *)datum {
    ReferenceStackEntryMO *stackEntryMO = [NSEntityDescription insertNewObjectForEntityForName:@"ReferenceStackEntry" inManagedObjectContext:context];
    stackEntryMO.relationship = @"reference" ;
    stackEntryMO.datum = datum ;
    stackEntryMO.sequence = sequence ;



    return stackEntryMO ;
}
@end

NSString *generateSnpString(NSSet *set) {
    NSString *datumSnpString = @"";
    for(DatumSnpMO *datumSnpMO in set){
        datumSnpString = [datumSnpString stringByAppendingFormat:@"%ld,",datumSnpMO.column.integerValue];
    }
    if(datumSnpString.length>0){
        datumSnpString = [datumSnpString substringToIndex:datumSnpString.length-1];
    }
    return datumSnpString;
}