//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//
//
//


#import "SampleRepository.h"
#import "SampleMO.h"


@implementation SampleRepository {

}


+ (SampleRepository *)sharedInstance
{
    static SampleRepository *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[SampleRepository alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}

- (SampleMO *)getSampleForName:(NSString *)name andContext:(NSManagedObjectContext *)moc andError:(NSError *)error {
    NSEntityDescription *entityDescription = [NSEntityDescription
            entityForName:@"Sample" inManagedObjectContext:moc];
    NSFetchRequest *request = [[NSFetchRequest alloc] init];
    [request setEntity:entityDescription];
    NSPredicate *predicate = [NSPredicate predicateWithFormat:@"name == %@", name];
    [request setPredicate:predicate];
    NSArray *sampleArray = [moc executeFetchRequest:request error:&error];
    SampleMO *sampleMO = nil ;
    if (error || sampleArray == nil || sampleArray.count == 0) {
        if(sampleArray==nil){
            NSLog(@"sampleArray is nil: %@ ", error);
        }
        else{
            NSLog(@"sampleArray %ld error %@", sampleArray.count, error);
        }
        return nil ;
    }
    else {
        sampleMO = [sampleArray objectAtIndex:0];
//                    NSLog(@"sample found %@",sampleMO.name);
    }
    return sampleMO;
}


- (SampleMO *)insertSample:(NSManagedObjectContext *)context id:(NSNumber *)number name:(NSString *)name{
    SampleMO *sampleMO = [NSEntityDescription insertNewObjectForEntityForName:@"Sample" inManagedObjectContext:context];
    sampleMO.sampleId = number ;
    sampleMO.name = name ;
    return sampleMO;
}

- (NSArray *)getAllSamples:(NSManagedObjectContext *)context {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Sample" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];
    NSError *error1;
    return [context executeFetchRequest:request1 error:&error1];
}

@end