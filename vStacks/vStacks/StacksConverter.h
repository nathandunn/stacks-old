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
//@class DepthRepository;
//@class HaplotypeRepository;
@class LocusRepository;
//@class SnpRepository;
//@class StackEntryRepository;
//@class AlleleRepository;
@class ProgressController;


@interface StacksConverter : NSObject


//@property(nonatomic, strong) DatumRepository* datumRepository ;
//@property(nonatomic, strong) DepthRepository* depthRepository;
//@property(nonatomic, strong) HaplotypeRepository* haplotypeRepository;
//@property(nonatomic, strong) LocusRepository* locusRepository;
//@property(nonatomic, strong) PopulationRepository* populationRepository;
//@property(nonatomic, strong) SampleRepository* sampleRepository ;
//@property(nonatomic, strong) SnpRepository* snpRepository ;
//@property(nonatomic, strong) StackEntryRepository* stackEntryRepository ;
//@property(nonatomic, strong) AlleleRepository* alleleRepository;

// a lookup
//@property(nonatomic, strong) NSMutableDictionary *lociDictionary ;
// sample:Dictionary<internalid,externalid>
//@property(nonatomic, strong) NSMutableDictionary *sampleLookupDictionary;
//@property(nonatomic, strong) NSMutableDictionary *locusSnpMap;
@property(nonatomic, strong) NSNumberFormatter *numberFormatter;


@property(nonatomic) bool stopProcess;

@property(atomic, strong) NSPersistentStoreCoordinator *persistentStoreCoordinator ;

- (id)init ;

- (NSString *) generateFilePathForUrl:(NSURL *) url;

- (StacksDocument *)loadLociAndGenotypes:(NSString *)path progressWindow:(ProgressController *)progressController importPath:(NSString *)path1;
- (NSMutableDictionary *)loadPopulation:(NSString *)path;

- (StacksDocument *)loadDocument:(StacksDocument *)document progressWindow:(ProgressController *)bar importPath:(NSString *)path;
- (StacksDocument *)createStacksDocumentForPath:(NSString *)path;

//- (NSManagedObjectContext *)getContextForPath:(NSString *)path andDocument:(StacksDocument *) document;
- (NSManagedObjectContext *)getContextForPath:(NSString *)string andName:(NSString *)name  andDocument:(StacksDocument *) document;


//- (StacksDocument *)getStacksDocumentForPath:(NSString *)string;
//- (void)loadLociAndGenotypes:(NSString *)appendingString progressBar:(NSProgressIndicator *)bar;

- (void)addPopulationsToDocument:(StacksDocument *)document forPath:(NSString *)path ;


@end

NSUInteger  countParents(NSArray *parents);


