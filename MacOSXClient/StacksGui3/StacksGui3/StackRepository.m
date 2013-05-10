//
// Created by Nathan Dunn on 5/8/13.
// Copyright (c) 2013 Nathan Dunn. All rights reserved.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StackRepository.h"
#import "StackMO.h"
#import "DatumMO.h"


@implementation StackRepository {

}
- (StackMO *)insertStack:(NSManagedObjectContext *)context datum:(DatumMO *)datum {
    StackMO *stackMO = [NSEntityDescription insertNewObjectForEntityForName:@"Stack" inManagedObjectContext:context];
    stackMO.datum = datum;
    datum.stack = stackMO;
    return stackMO;
}

- (StackMO *)getStack:(NSManagedObjectContext *)context forDatum:(DatumMO *)datum {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Stack" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];

    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"datum == %@ ", datum];
    [request1 setPredicate:predicate1];
    NSError *error1;
    NSArray *locusArray = [context executeFetchRequest:request1 error:&error1];
    if(locusArray!=nil && locusArray.count==1){
        return [locusArray objectAtIndex:0] ;
    }
    else{
        return nil ;
    }
}
@end