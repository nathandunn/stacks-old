//
// Created by NathanDunn on 2/28/13.
//
//
//


#import <Foundation/Foundation.h>

@class StacksDocument;
@class SampleRepository;
@class PopulationRepository;
@class DatumRepository;
@class LocusRepository;
@class ProgressController;


@interface StacksConverter : NSObject


@property(nonatomic, strong) NSNumberFormatter *numberFormatter;


@property(nonatomic) bool stopProcess;

@property(atomic, strong) NSPersistentStoreCoordinator *persistentStoreCoordinator;

- (id)init;

- (NSString *)generateFilePathForUrl:(NSURL *)url;

- (StacksDocument *)loadLociAndGenotypes:(NSString *)path progressWindow:(ProgressController *)progressController importPath:(NSString *)path1;

- (NSMutableDictionary *)loadPopulation:(NSString *)path;

- (StacksDocument *)loadDocument:(StacksDocument *)document progressWindow:(ProgressController *)bar importPath:(NSString *)path;

- (StacksDocument *)createStacksDocumentForPath:(NSString *)path;

- (NSManagedObjectContext *)getContextForPath:(NSString *)string andName:(NSString *)name andDocument:(StacksDocument *)document;

- (void)addPopulationsToDocument:(StacksDocument *)document forPath:(NSString *)path;


@end

NSUInteger countParents(NSArray *parents);


