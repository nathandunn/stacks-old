//
//  main.m
//  StacksGui3
//
//  Created by Nathan Dunn on 4/18/13.
//  Copyright (c) 2013 Nathan Dunn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

//NSManagedObjectContext *managedObjectContext();
//NSManagedObjectModel *managedObjectModel();
//NSURL *applicationLogDirectory() ;

int main(int argc, char *argv[])
{
    return NSApplicationMain(argc, (const char **)argv);
}

/*
NSManagedObjectContext *managedObjectContext()
{
    static NSManagedObjectContext *moc = nil;
    
    if (moc != nil) {
        return moc;
    }
    
    NSPersistentStoreCoordinator *coordinator =
    [[NSPersistentStoreCoordinator alloc]
     initWithManagedObjectModel: managedObjectModel()];
    
    NSString *STORE_TYPE = NSXMLStoreType;
    NSString *STORE_FILENAME = @"CDCLI.cdcli";
    
    NSError *error;
    NSURL *url = [applicationLogDirectory() URLByAppendingPathComponent:STORE_FILENAME];
    
    NSPersistentStore *newStore = [coordinator addPersistentStoreWithType:STORE_TYPE
                                                            configuration:nil URL:url options:nil
                                                                    error:&error];
    
    if (newStore == nil) {
        
        NSLog(@"Store Configuration Failure\n%@",
              ([error localizedDescription] != nil) ?
              [error localizedDescription] : @"Unknown Error");
    }
    
    
    moc = [[NSManagedObjectContext alloc] initWithConcurrencyType:NSMainQueueConcurrencyType];
    [moc setPersistentStoreCoordinator:coordinator];
    
    return moc;
}

NSManagedObjectModel *managedObjectModel() {
    
    static NSManagedObjectModel *mom = nil;
    
    if (mom != nil) {
        return mom;
    }
    
    NSEntityDescription *runEntity = [[NSEntityDescription alloc] init];
    [runEntity setName:@"Run"];
    [runEntity setManagedObjectClassName:@"Run"];
    
    NSAttributeDescription *dateAttribute = [[NSAttributeDescription alloc] init];
    
    [dateAttribute setName:@"date"];
    [dateAttribute setAttributeType:NSDateAttributeType];
    [dateAttribute setOptional:NO];
    
    
    NSAttributeDescription *idAttribute = [[NSAttributeDescription alloc] init];
    
    [idAttribute setName:@"processID"];
    [idAttribute setAttributeType:NSInteger64AttributeType];
    [idAttribute setOptional:NO];
    [idAttribute setDefaultValue:@(-1)];
    
    NSExpression *lhs = [NSExpression expressionForEvaluatedObject];
    NSExpression *rhs = [NSExpression expressionForConstantValue:@0];
    
    NSPredicate *validationPredicate = [NSComparisonPredicate
                                        predicateWithLeftExpression:lhs
                                        rightExpression:rhs
                                        modifier:NSDirectPredicateModifier
                                        type:NSGreaterThanPredicateOperatorType
                                        options:0];
    
    NSString *validationWarning = @"Process ID < 1";
    
    [idAttribute setValidationPredicates:@[validationPredicate]
                  withValidationWarnings:@[validationWarning]];
    
    [runEntity setProperties:@[dateAttribute, idAttribute]];
    
    mom = [[NSManagedObjectModel alloc] init];
    [mom setEntities:@[runEntity]];
    
    NSDictionary *localizationDictionary = @{
                                             @"Property/date/Entity/Run":@"Date",
                                             @"Property/processID/Entity/Run":@"Process ID",
                                             @"ErrorString/Process ID < 1":@"Process ID must not be less than 1"};
    
    [mom setLocalizationDictionary:localizationDictionary];
    
    return mom;
}

NSURL *applicationLogDirectory() {
    
    NSString *LOG_DIRECTORY = @"CDCLI";
    static NSURL *ald = nil;
    
    if (ald == nil) {
        
        NSFileManager *fileManager = [[NSFileManager alloc] init];
        NSError *error;
        NSURL *libraryURL = [fileManager URLForDirectory:NSLibraryDirectory inDomain:NSUserDomainMask appropriateForURL:nil create:YES error:&error];
        if (libraryURL == nil) {
            NSLog(@"Could not access Library directory\n%@", [error localizedDescription]);
        }
        else {
            ald = [libraryURL URLByAppendingPathComponent:@"Logs"];
            ald = [ald URLByAppendingPathComponent:LOG_DIRECTORY];
            NSDictionary *properties = [ald resourceValuesForKeys:@[NSURLIsDirectoryKey]
                                                            error:&error];
            if (properties == nil) {
                if (![fileManager createDirectoryAtURL:ald withIntermediateDirectories:YES attributes:nil error:&error]) {
                    NSLog(@"Could not create directory %@\n%@", [ald path], [error localizedDescription]);
                    ald = nil;
                }
            }
        }
    }
    return ald;
}
*/