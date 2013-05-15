//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class StacksView;
@class LocusView;
@class StacksDocument;
@class SampleRepository;
@class PopulationRepository;
@class DatumRepository;
@class DepthRepository;
@class HaplotypeRepository;
@class LocusRepository;
@class SnpRepository;
@class StackEntryRepository;


@interface StacksConverter : NSObject


@property(nonatomic, strong) DatumRepository* datumRepository ;
@property(nonatomic, strong) DepthRepository* depthRepository;
@property(nonatomic, strong) HaplotypeRepository* haplotypeRepository;
@property(nonatomic, strong) LocusRepository* locusRepository;
@property(nonatomic, strong) PopulationRepository* populationRepository;
@property(nonatomic, strong) SampleRepository* sampleRepository ;
@property(nonatomic, strong) SnpRepository* snpRepository ;
@property(nonatomic, strong) StackEntryRepository* stackEntryRepository ;
//@property(nonatomic, strong) NSManagedObjectContext *managedObjectContext;
//@property(nonatomic, strong) NSManagedObjectModel *managedObjectModel;

- (id)init ;
- (StacksView *)loadStacksView:(NSString *)filename atPath:(NSString *)path forTag:(NSInteger)tag locus:(LocusView *)locus;
- (StacksDocument *)loadLociAndGenotypes:(NSString *)path;
- (NSMutableDictionary *)loadPopulation:(NSString *)path;

- (StacksDocument *)loadDocument:(StacksDocument *)document;

- (StacksDocument *)createStacksDocumentForPath:(NSString *)path;

- (StacksDocument *)getStacksDocumentForPath:(NSString *)string;
@end


