//
// Created by Nathan Dunn on 4/1/14.
// Copyright (c) 2014 Nathan Dunn. All rights reserved.
//

#import "GenericHashRepository.h"
#import "GenericHashMO.h"


@implementation GenericHashRepository {

}

+ (GenericHashRepository *)sharedInstance
{
    static GenericHashRepository *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[GenericHashRepository alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}

- (void)storeKey:(NSString *)sampleString dataValue:(NSMutableDictionary *)value type:(NSString *)type {

}

- (GenericHashMO *)store:(NSManagedObjectContext *)context key:(NSString *)key dictionary:(NSMutableDictionary *)dictionary type:(NSString *)type {
    NSData *myData = [NSKeyedArchiver archivedDataWithRootObject:dictionary];

    GenericHashMO *genericHasHMO = [NSEntityDescription insertNewObjectForEntityForName:@"GenericHash" inManagedObjectContext:context];
    genericHasHMO.key = key ;
    genericHasHMO.dataValue = myData ;
    genericHasHMO.type = type ;
    genericHasHMO.stringValue = nil ;

    return genericHasHMO;
}

- (NSDictionary *)getDictionary:(NSManagedObjectContext *)context forKey:(NSString *)key {
    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"GenericHash" inManagedObjectContext:context];
    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
    [request1 setEntity:entityDescription1];

    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"key == %@ ", key];
    [request1 setPredicate:predicate1];
    NSError *error1;
    NSArray *resultArray = [context executeFetchRequest:request1 error:&error1];
    if(resultArray !=nil && resultArray.count==1){
//        NSString* className = [[resultArray objectAtIndex:0] className];
        // is GenericHasMO
        GenericHashMO* resultData = [resultArray objectAtIndex:0] ;

        NSDictionary *myDictionary = (NSDictionary*) [NSKeyedUnarchiver unarchiveObjectWithData:resultData.dataValue];
        return myDictionary ;
    }
    else{
        return nil ;
    }
}

/**
* TODO: remove all
*/
- (void)detachAll {

}
@end