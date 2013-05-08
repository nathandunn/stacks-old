//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "SampleRepository.h"
#import "SampleMO.h"


@implementation SampleRepository {

}
- (SampleMO *)getSampleForName:(NSString *)name andContext:(NSManagedObjectContext *)moc andError:(NSError *)error {
//    if (error == nil) {
//        error = [[NSError alloc] init];
//    }
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


- (SampleMO *)insertSample:(NSManagedObjectContext *)context withId:(NSNumber *) number andName:(NSString *)name{
    SampleMO *sampleMO = [NSEntityDescription insertNewObjectForEntityForName:@"Sample" inManagedObjectContext:context];
    sampleMO.sampleId = number ;
    sampleMO.name = name ;
    return sampleMO;
}
@end