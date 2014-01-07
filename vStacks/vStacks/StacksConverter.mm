//
// Created by NathanDunn on 2/28/13.
//
//



#include "LocusMO.h"
//#include "stacks.h"
#include "sql_utilities.h"
#include "PopMap.h"
#import "StacksDocument.h"


#include <fstream>

using std::ifstream;
using std::ofstream;


//#import "StackEntryMO.h"
#import "SampleRepository.h"
#import "StacksConverter.h"
#import "PopulationRepository.h"
//#import "DepthRepository.h"
#import "DatumRepository.h"
//#import "HaplotypeRepository.h"
#import "LocusRepository.h"
#import "DatumMO.h"
//#import "HaplotypeMO.h"
//#import "DepthMO.h"
#import "PopulationMO.h"
#import "SampleMO.h"
//#import "SnpRepository.h"
//#import "StackEntryRepository.h"
//#import "DatumSnpMO.h"
//#import "DatumAlleleMO.h"
//#import "AlleleRepository.h"
//#import "LocusAlleleMO.h"
//#import "ConsensusStackEntryMO.h"
#import "ProgressController.h"
#import "StackEntryRenderer.h"


#include <sys/time.h>


#define CHECK_STOP  if (stopProcess) { [progressWindow close]; return nil ; }


NSString *calculateType(NSString *file);

//void setParentCounts(NSMutableDictionary *dictionary, NSString *file);
//void setParentCounts(NSString *file);


@implementation StacksConverter {

}


/**
* https://casspr.fogbugz.com/default.asp?1144
* 1 - add pops
2 - add samples
3 - load locis (refactor into other method)
4 - create datum from loci and sample (do per sample, load loci into a lookupTable to use throughout)
5 - add stackEntries onto each datum per tags file
6 - load snps onto each datum
7 - load alleles onto each datum
*/

// repositories
//@synthesize [DatumRepository sharedInstance];
//@synthesize [DepthRepository sharedInstance];
//@synthesize [HaplotypeRepository sharedInstance];
//@synthesize [LocusRepository sharedInstance];
//@synthesize [PopulationRepository sharedInstance];
//@synthesize [SampleRepository sharedInstance];
//@synthesize [SnpRepository sharedInstance];
//@synthesize stackEntryRepository;
//@synthesize [AlleleRepository sharedInstance];


// lookups
//@synthesize lociDictionary;
@synthesize sampleLookupDictionary;
@synthesize stopProcess;
@synthesize locusSnpMap;


@synthesize persistentStoreCoordinator;


- (id)init {
    self = [super init];
    if (self) {
        // nothing write now
//        [DatumRepository sharedInstance] = [[DatumRepository alloc] init];
//        [DepthRepository sharedInstance] = [[DepthRepository alloc] init];
//        [HaplotypeRepository sharedInstance] = [[HaplotypeRepository alloc] init];
//        [LocusRepository sharedInstance] = [[LocusRepository alloc] init];
//        [PopulationRepository sharedInstance] = [[PopulationRepository alloc] init];
//        [SampleRepository sharedInstance] = [[SampleRepository alloc] init];
//        [SnpRepository sharedInstance] = [[SnpRepository alloc] init];
//        stackEntryRepository = [[StackEntryRepository alloc] init];
//        [AlleleRepository sharedInstance] = [[AlleleRepository alloc] init];


//        lociDictionary = [[NSMutableDictionary alloc] init];
        sampleLookupDictionary = [[NSMutableDictionary alloc] init];
        stopProcess = false;


        persistentStoreCoordinator = [[NSPersistentStoreCoordinator alloc] initWithManagedObjectModel:[NSManagedObjectModel mergedModelFromBundles:nil]];
    }
    return self;
}

- (NSString *)generateFilePathForUrl:(NSURL *)url {
    return [url.path stringByAppendingFormat:@"/%@.stacks", url.path.lastPathComponent];
}


- (StacksDocument *)loadLociAndGenotypes:(NSString *)path progressWindow:(ProgressController *)progressController {

    StacksDocument *stacksDocument = [self createStacksDocumentForPath:path];
    if (stacksDocument == nil) {
        return nil;
    }

    stacksDocument = [self loadDocument:stacksDocument progressWindow:progressController];
    NSManagedObjectContext *moc = stacksDocument.managedObjectContext;
    NSError *error;
    [moc save:&error];
    NSLog(@"error %@", error);
    [[moc parentContext] save:&error];
    NSLog(@"error 2 %@", error);
    [moc reset];
//    [[moc parentContext] save:&error];
//    [[moc parentContext] performBlock:^(){
//        [[moc parentContext] save:NULL];
//    }];

//    [lociDictionary removeAllObjects];
    [sampleLookupDictionary removeAllObjects];
    [moc reset];
    [[moc parentContext] reset];


    return stacksDocument;
}

//- (NSManagedObjectContext *)getContextForPath:(NSString *)path andDocument:(StacksDocument *)document{
//    if(document.name==NULL){
//        return [self getContextForPath:path andName:@"StacksDocument" andDocument:document];
//    }
//    else{
//
//        return [self getContextForPath:path andName:document.name  andDocument:document];
//    }
//}


- (NSManagedObjectContext *)getContextForPath:(NSString *)path andName:(NSString *)name andDocument:(StacksDocument *)document {
    NSDictionary *options = [NSDictionary dictionaryWithObjectsAndKeys:[NSNumber numberWithBool:YES], NSMigratePersistentStoresAutomaticallyOption,
                                                                       [NSNumber numberWithBool:YES],
                                                                       NSInferMappingModelAutomaticallyOption, nil];


    NSURL *storeUrl = [NSURL fileURLWithPath:[path stringByAppendingFormat:@"/%@.stacks", name]];
    NSLog(@"saving to %@ from %@", path, storeUrl);
    if (persistentStoreCoordinator == nil) {
        persistentStoreCoordinator = [[NSPersistentStoreCoordinator alloc] initWithManagedObjectModel:[NSManagedObjectModel mergedModelFromBundles:nil]];
    }
    NSError *error = nil;

    if (![persistentStoreCoordinator addPersistentStoreWithType:NSSQLiteStoreType configuration:nil URL:storeUrl options:options error:&error]) {
        NSLog(@"error loading persistent store..");
        [[NSFileManager defaultManager] removeItemAtPath:storeUrl.path error:nil];
        if (![persistentStoreCoordinator addPersistentStoreWithType:NSSQLiteStoreType configuration:nil URL:storeUrl options:options error:&error]) {
            NSLog(@"Unresolved error %@, %@", error, [error userInfo]);
            //abort();
        }
    }


    NSManagedObjectContext *context = [[NSManagedObjectContext alloc] init];
    [context setPersistentStoreCoordinator:persistentStoreCoordinator];
    document.managedObjectContext = context;
    document.path = path;


    return context;
}


- (NSMutableDictionary *)loadPopulation:(NSString *)path {
    NSMutableDictionary *populationLookup = [[NSMutableDictionary alloc] init];

    NSString *popmapFile = [path stringByAppendingString:@"popmap"];
    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL exists = [fileManager fileExistsAtPath:popmapFile];

    if (exists) {
        NSArray *fileData = [[NSString stringWithContentsOfFile:popmapFile encoding:NSUTF8StringEncoding error:nil] componentsSeparatedByString:@"\n"];
        NSString *line;
        for (line in fileData) {
            NSArray *columns = [line componentsSeparatedByString:@"\t"];
            if (columns.count == 2) {
                [populationLookup setObject:[columns objectAtIndex:1] forKey:[columns objectAtIndex:0]];
            }
            else {
                NSLog(@"something wrong with the column count %@", line);
            }
        }
    }
    else {
        NSLog(@"does not exist at %@", popmapFile);
    }

    return populationLookup;
}

/**
* returns all of the files that end with .tags.tsv in the directory
*/
- (vector<pair<int, string> >)buildFileList:(NSString *)path {
    vector<pair<int, string> > files;

    NSFileManager *fileManager = [NSFileManager defaultManager];

    NSDirectoryEnumerator *dirEnum = [fileManager enumeratorAtPath:path];
    NSString *file;
    int i = 0;
    while (file = [dirEnum nextObject]) {

        if ([file hasSuffix:@".tags.tsv"]) {
//            NSLog(@"HAS SUFFIX name %@", file);
            NSUInteger length = [file length] - 9;
            NSString *fileName = [file substringToIndex:length];
            files.push_back(make_pair(i, [fileName UTF8String]));
            ++i;
        }
    }

    sort(files.begin(), files.end());

    return files;
}

/**
* Exist on error.
*/
- (NSString *)checkFile:(NSString *)examplePath {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];

    if (!existsAtPath) {
        NSLog(@"files do not exist %@", examplePath);
        exit(1);
    }

    NSError *error;
    NSArray *files = [fileManager contentsOfDirectoryAtPath:examplePath error:&error];
    if (error) {
        NSLog(@"error %@", error);
        return nil;
    }
    NSPredicate *fltr = [NSPredicate predicateWithFormat:@"self ENDSWITH '.catalog.tags.tsv'"];
    NSPredicate *batch = [NSPredicate predicateWithFormat:@"self BEGINSWITH 'batch'"];
    NSArray *onlyCatalog = [[files filteredArrayUsingPredicate:fltr] filteredArrayUsingPredicate:batch];

    NSLog(@"# of files: %ld", onlyCatalog.count);
    for (NSString *file in onlyCatalog) {
        NSLog(@"file %@", file);
    }

    NSString *originalFile = [onlyCatalog objectAtIndex:0];
    NSUInteger firstIndex = [originalFile rangeOfString:@"."].location;
    NSString *batchName = [originalFile substringToIndex:firstIndex];

    return batchName;

}

- (StacksDocument *)loadDocument:(StacksDocument *)stacksDocument progressWindow:(ProgressController *)progressWindow {
    NSProgressIndicator *bar = progressWindow.loadProgress;
    if (bar != nil) {
        progressWindow.actionTitle.stringValue = [NSString stringWithFormat:@"Loading %@", stacksDocument.name];
        progressWindow.actionMessage.stringValue = @"Begin import";
        bar.doubleValue = 0;
        [bar display];
//        [bar incrementBy:1];
    }
    NSString *path = stacksDocument.path;
    // returns batch_1, batch_2, etc. whatever exists
    NSString *batchName = [self checkFile:path];
    NSLog(@"returned batch name: %@", batchName);
    map<int, CSLocus *> catalog;
    NSString *catalogTagFile = [path stringByAppendingFormat:@"%@.catalog.tags.tsv", batchName];
    stacksDocument.type = calculateType(catalogTagFile);
    NSLog(@"type is %@", stacksDocument.type);
    NSString *catalogFile = [path stringByAppendingFormat:@"%@.catalog", batchName];

    struct timeval time1, time2;
    gettimeofday(&time1, NULL);


    load_loci([catalogFile UTF8String], catalog, false);
    progressWindow.actionMessage.stringValue = @"Loading loci";
//    [bar incrementBy:2];
    gettimeofday(&time2, NULL);
    NSLog(@"load_loci %ld", (time2.tv_sec - time1.tv_sec));

    /**
    * START: loading datums + extra info
*/
    vector<vector<CatMatch *> > catalog_matches;
    map<int, string> samples;
    vector<int> sample_ids;

    progressWindow.actionMessage.stringValue = @"Building file list";
    vector<pair<int, string>> files = [self buildFileList:path];
//    [bar incrementBy:2];
    NSLog(@"number of files %ld", files.size());



    // loci loaded . . . now loading datum
    progressWindow.actionMessage.stringValue = @"Matching catalogs";
    NSLog(@"model size %d", (int) catalog.size());

    double incrementAmount = 5.0 / files.size();
    gettimeofday(&time1, NULL);
    for (uint i = 0; i < files.size(); i++) {
        vector<CatMatch *> m;
        NSString *sampleString = [NSString stringWithUTF8String:files[i].second.c_str()];
        NSString *matchString = [path stringByAppendingFormat:@"/%@", sampleString];
//        NSLog(@"loading match file %@ for sample name %@", matchString, sampleString);
        if (([sampleString rangeOfString:@"catalog"]).location == NSNotFound) {
            NSMutableDictionary *matchDictionary = [self loadMatchesDictionary:[matchString stringByAppendingString:@".matches.tsv"]];
            [sampleLookupDictionary setObject:matchDictionary forKey:sampleString];
        }

        load_catalog_matches([matchString UTF8String], m);

        if (m.size() == 0) {
            cerr << "Warning: unable to find any matches in file '" << files[i].second << "', excluding this sample from population analysis.\n";
            continue;
        }

        catalog_matches.push_back(m);
        if (samples.count(m[0]->sample_id) == 0) {
            samples[m[0]->sample_id] = files[i].second;
            sample_ids.push_back(m[0]->sample_id);
        } else {
            cerr << "Fatal error: sample ID " << m[0]->sample_id << " occurs twice in this stacksDocuments set, likely the pipeline was run incorrectly.\n";
            exit(1);
        }
        [bar incrementBy:incrementAmount];
    }
    gettimeofday(&time2, NULL);
    NSLog(@"catalog matches %ld", (time2.tv_sec - time1.tv_sec));


    NSLog(@"input populations %ld", stacksDocument.populations.count);
    progressWindow.actionMessage.stringValue = @"Adding populations";
    [self addPopulationsToDocument:stacksDocument forPath:path];
    [bar incrementBy:5];
    progressWindow.actionMessage.stringValue = @"Adding samples";
    NSLog(@"output populations %ld", stacksDocument.populations.count);
    [self addSamplesToDocument:stacksDocument forSampleIds:sample_ids andSamples:samples];
    [bar incrementBy:5];
    NSLog(@"output populations after samples %ld", stacksDocument.populations.count);
    for (PopulationMO *populationMo in stacksDocument.populations) {
        NSLog(@"samples %ld per population %@", populationMo.samples.count, populationMo.name);
//        for(SampleMO* sampleMO in populationMo.samples){
//            NSLog(@"sample %@ ",sampleMO.name);
//        }
    }

    cerr << "Populating observed haplotypes for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";

    gettimeofday(&time1, NULL);
    progressWindow.actionMessage.stringValue = @"Populating samples";
    PopMap<CSLocus> *pmap = new PopMap<CSLocus>((int) sample_ids.size(), (int) catalog.size());
    pmap->populate(sample_ids, catalog, catalog_matches);

//    delete catalog_matches;
//    catalog_matches.erase (catalog_matches.begin(),catalog_matches.end());


    [bar incrementBy:5];
    gettimeofday(&time2, NULL);
    NSLog(@"population pmap %ld", (time2.tv_sec - time1.tv_sec));
    progressWindow.actionMessage.stringValue = @"Populated Samples";


    NSMutableSet *loci = [[NSMutableSet alloc] init];

    gettimeofday(&time1, NULL);
    map<int, CSLocus *>::iterator catalogIterator = catalog.begin();
    NSManagedObjectContext *moc = stacksDocument.managedObjectContext;

    NSNumberFormatter *numberFormatter = [[NSNumberFormatter alloc] init];
    numberFormatter.numberStyle = NSNumberFormatterNoStyle;

    incrementAmount = 10.0 / catalog.size();

    progressWindow.actionMessage.stringValue = @"Loading locus snps";
    while (catalogIterator != catalog.end()) {
        const char *read = (*catalogIterator).second->con;
        LocusMO *locusMO = [[LocusRepository sharedInstance] insertNewLocus:moc withId:[NSNumber numberWithInt:(*catalogIterator).second->id]
                                                               andConsensus:[[NSString alloc] initWithCString:read encoding:NSUTF8StringEncoding] andMarker:[NSString stringWithUTF8String:catalogIterator->second->marker.c_str()]
        ];

        // get catalogs for matches
        // TODO: double-check that this is correct . . .
        locusMO.length = [NSNumber numberWithInt:catalogIterator->second->depth];
        locusMO.type = stacksDocument.type;

//        NSLog(@"chromosme %@",[NSString stringWithUTF8String:catalogIterator->second->loc.chr]);

        locusMO.chromosome = [NSString stringWithUTF8String:catalogIterator->second->loc.chr];
        unsigned int intValue = (unsigned int) catalogIterator->second->loc.bp;
//        NSLog(@"int value %@",[NSNumber numberWithUnsignedInt:intValue]);
//        locusMO.basePairs = [NSNumber numberWithUnsignedInt:catalogIterator->second->loc.bp];
        locusMO.basePairs = [NSNumber numberWithUnsignedInt:intValue];
        locusMO.strand = catalogIterator->second->loc.strand == plus ? @"+" : @"-";

//        catalogIterator->second->
        vector<SNP *> snps = catalogIterator->second->snps;
//        NSLog(@"inserting locus %@ with sequence %@", locusMO.locusId, [NSString stringWithUTF8String:read]);

        vector<SNP *>::iterator snpsIterator = snps.begin();


        NSMutableArray *snpArray = [[NSMutableArray alloc] init];
        for (; snpsIterator != snps.end(); ++snpsIterator) {
            SNP *snp = (*snpsIterator);

            NSDictionary *snpDictionary = [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSString stringWithFormat:@"%i", snp->col], @"column"
                    , [NSNumber numberWithFloat:snp->lratio], @"lratio"
                    , [NSNumber numberWithFloat:snp->rank_1], @"rank1"
                    , [NSNumber numberWithFloat:snp->rank_2], @"rank2"
                    , [NSNumber numberWithFloat:snp->rank_3], @"rank3"
                    , [NSNumber numberWithFloat:snp->rank_4], @"rank4"
                    , nil ];
            [snpArray addObject:snpDictionary];

//            NSError* error ;
//            NSData* snpData = [NSJSONSerialization dataWithJSONObject:snpDictionary options:NSJSONWritingPrettyPrinted error:&error];
//            NSString* snpString = [[NSString alloc] initWithData:snpData encoding:NSUTF8StringEncoding];
//
//            NSLog(@"snpString %@",snpString) ;

        }

        NSError *error2;
        NSData *snpArrayData = [NSJSONSerialization dataWithJSONObject:snpArray options:NSJSONWritingPrettyPrinted error:&error2];
        locusMO.snpData = snpArrayData;
//        NSString *snpArrayString = [[NSString alloc] initWithData:snpArrayData encoding:NSUTF8StringEncoding];
//        NSLog(@"snpArrayString %@",snpArrayString) ;

//        map<string, int> alleles;   // Map of the allelic configuration of SNPs in this stack along with the count of each
        map<string, int> alleles = catalogIterator->second->alleles;
        map<string, int>::iterator allelesIterator = alleles.begin();
//        string allele;
//        int column;
        NSMutableArray *alleleArray = [[NSMutableArray alloc] init];
        for (; allelesIterator != alleles.end(); ++allelesIterator) {
            string allele = allelesIterator->first;
            int column = allelesIterator->second;

            NSDictionary *alleleDictionary = [NSDictionary dictionaryWithObjectsAndKeys:
                    [NSNumber numberWithInt:column], @"depth"
                    , [NSString stringWithUTF8String:allele.c_str()], @"allele"
//                    ,[NSNumber numberWithFloat:snp->rank_1],@"ratio"
                    , nil ];
            [alleleArray addObject:alleleDictionary];
        }

//        NSError *error3;
        NSData *alleleArrayData = [NSJSONSerialization dataWithJSONObject:alleleArray options:NSJSONWritingPrettyPrinted error:&error2];
//        NSString *alleleArrayString = [[NSString alloc] initWithData:alleleArrayData encoding:NSUTF8StringEncoding];
//        NSLog(@"alleleArrayString %@",alleleArrayString) ;
        locusMO.alleleData = alleleArrayData;


        [loci addObject:locusMO];
//        [lociDictionary setObject:locusMO forKey:locusMO.locusId];
        ++catalogIterator;
        [bar incrementBy:incrementAmount];
    }
    gettimeofday(&time2, NULL);


//    setParentCounts(lociDictionary, catalogTagFile);
    [self setParentCounts:stacksDocument.managedObjectContext forFile:catalogTagFile];

    NSLog(@"populating snps %ld", (time2.tv_sec - time1.tv_sec));

    map<int, CSLocus *>::iterator it;
    Datum *datum;
    CSLocus *loc;


    CHECK_STOP
    progressWindow.actionMessage.stringValue = @"Loading datums";
    // for each sample process the catalog
    NSLog(@"samples %ld X catalog %ld = %ld ", sample_ids.size(), catalog.size(), sample_ids.size() * catalog.size());

    long totalCatalogTime = 0;
    incrementAmount = 30.0 / (sample_ids.size() * catalog.size());

    // 7 is 400 X 7 = 3K . . .
    uint saveAfterSamples = 100000;


    gettimeofday(&time1, NULL);

    long iterCount = 0;
    //go through all samples
    for (uint i = 0; i < sample_ids.size(); i++) {
        int sampleId = sample_ids[i];
        CHECK_STOP
        string sampleString = samples[sampleId];
        progressWindow.actionMessage.stringValue = [NSString stringWithFormat:@"Loading datum %i/%ld", i + 1, sample_ids.size()];

        NSString *key = [NSString stringWithUTF8String:sampleString.c_str()];
//        SampleMO *sampleMO = [[SampleRepository sharedInstance] getSampleForName:key andContext:moc andError:nil];

        gettimeofday(&time1, NULL);
        // go through all loci
        for (it = catalog.begin(); it != catalog.end(); it++, iterCount++) {
            loc = it->second;
            datum = pmap->datum(loc->id, sample_ids[i]);

            // TODO: use a lookup here to speed up

//            NSNumber *lookupKey = [NSNumber numberWithInt:it->first];
//            NSNumber *lookupKey = [NSNumber numberWithInteger:[[NSString stringWithFormat:@"%d", it->first] integerValue]];
//            locusMO = [lociDictionary objectForKey:[NSString stringWithFormat:@"%ld", it->first]];
//            LocusMO *locusMO = [[LocusRepository sharedInstance] getLocus:stacksDocument.managedObjectContext forId:lookupKey.integerValue];
//            locusMO = [lociDictionary objectForKey:lookupKey];
//            if (locusMO == nil) {
//                for (int i = 0; i < 10; i++) {
//                    NSLog(@"Could not find Locus! %d", it->first);
//                }
//                return nil;
//            }

            if (datum != NULL) {

                vector<char *> obshape = datum->obshap;
                vector<int> depths = datum->depth;
                int numLetters = (int) obshape.size();
                // TODO: should be entering this for the locus as well?
//                if (loc->id == 1) {
//                    NSLog(@"insertign datum for sample %@ locus %@ and sampleId %i", sampleMO.name, locusMO.locusId, sample_ids[i]);
//                }
//                DatumMO *newDatumMO = [[DatumRepository sharedInstance] insertDatum:moc name:key sampleId:[NSNumber numberWithInt:sample_ids[i]] sample:sampleMO locus:locusMO];

                // TODO: NECESSARY, or can be done above??
                // TODO: is not the locus the same thing as the id?  can I use the loc->id here?
                DatumMO *newDatumMO = [NSEntityDescription insertNewObjectForEntityForName:@"Datum" inManagedObjectContext:moc];
                newDatumMO.name = key;
                newDatumMO.sampleId = [NSNumber numberWithInt:sampleId];
//                newDatumMO.sample = sample ;
//                newDatumMO.locus = locus ;
//                newDatumMO.tagId;
//                return newDatumMO ;
                newDatumMO.tagId = [NSNumber numberWithInt:datum->id];

                if (newDatumMO.sampleId == nil) {
                    NSLog(@"loading sample ID %@", newDatumMO.sampleId);
                }

                // get catalogs for matches
//                locusMO.length = [NSNumber numberWithInt:loc->depth];

                // TODO: CONVERT TO USE DATA
                NSMutableArray *datumDataArray = [[NSMutableArray alloc] init];
                if (depths.size() == numLetters && numLetters > 0 ) {
                    for (int j = 0; j < numLetters; j++) {

                        NSDictionary *dataDictionary = [NSDictionary dictionaryWithObjectsAndKeys:
                                [NSString stringWithUTF8String:obshape[j]], @"haplotype"
                                , [NSNumber numberWithInt:j], @"order"
                                , [NSNumber numberWithInt:depths[j]], @"depth"
                                , nil ];
//                        NSLog(@"haplotype %@",dataDictionary);
                        [datumDataArray addObject:dataDictionary];
//                        HaplotypeMO *haplotypeMO = [[HaplotypeRepository sharedInstance] insertHaplotype:moc haplotype:[NSString stringWithUTF8String:obshape[j]] andOrder:j];
//                        [newDatumMO addHaplotypesObject:haplotypeMO];
//
//                        DepthMO *depthMO = [[DepthRepository sharedInstance] insertDepth:moc depth:[NSNumber numberWithInt:depths[j]] andOrder:j];
//                        [newDatumMO addDepthsObject:depthMO];
                    }

                    NSError *error ;
                    NSData *datumData = [NSJSONSerialization dataWithJSONObject:datumDataArray options:NSJSONWritingPrettyPrinted error:&error];
                    newDatumMO.haplotypeData = datumData ;
//                    [locusMO addDatumsObject:newDatumMO];
                }
                else {
                    NSLog(@"mismatchon %@", [NSString stringWithUTF8String:sampleString.c_str()]);
                }
            }

            // end of process loci from catalogs

            // TODO: save within here!
            if (iterCount % saveAfterSamples == 0) {
                NSError *innerError = nil ;
                NSLog(@"saving samples");
                [stacksDocument.managedObjectContext save:&innerError];
                if (innerError != nil) {
                    NSLog(@"error doing inner save: %@", innerError);
                    return nil;
                }
            }
        }

        gettimeofday(&time2, NULL);
//        NSLog(@"iterating sample %d - time %ld", sample_ids[i], (time2.tv_sec - time1.tv_sec));
        totalCatalogTime += time2.tv_sec - time1.tv_sec;



//        else{
//            NSLog(@"NOT saving sample %i vs %i",i, (i%saveAfterSamples) );
//        }

        [bar incrementBy:incrementAmount];

    }

    NSLog(@"datum count: %ld", [[DatumRepository sharedInstance] getAllDatum:moc].count);


    NSError *innerError = nil ;
    stacksDocument.loci = loci;
    CHECK_STOP
    progressWindow.actionMessage.stringValue = @"Saving doc";
    [stacksDocument.managedObjectContext save:&innerError];
    [bar incrementBy:5];
    if (innerError != nil) {
        NSLog(@"error doing inner save: %@", innerError);
        return nil;
    }

    NSLog(@"total time %ld", totalCatalogTime);

    gettimeofday(&time1, NULL);
    delete pmap;
    gettimeofday(&time2, NULL);
    NSLog(@"delete pmap time %ld", time2.tv_sec - time1.tv_sec);


    gettimeofday(&time1, NULL);
    // TODO: I don't think this does anything hear
    CHECK_STOP
    progressWindow.actionMessage.stringValue = @"Reading populations";
    [self readPopulations:stacksDocument];

//    LocusMO *bLocusMO = [loci.allObjects objectAtIndex:0];
//    NSLog(@"pre locus %@ datums %ld", bLocusMO.locusId, bLocusMO.datums.count);
//
    stacksDocument.loci = loci;
//    LocusMO *cLocusMO = [stacksDocument.loci.allObjects objectAtIndex:0];
//    NSLog(@"post locus %@ datums %ld", cLocusMO.locusId, cLocusMO.datums.count);

    gettimeofday(&time2, NULL);
    NSLog(@"create stacks document time %ld", time2.tv_sec - time1.tv_sec);


    NSLog(@"loading stack entries");
    CHECK_STOP

//    for (id key in [sampleLookupDictionary allKeys]) {
//        NSLog(@"keys 1: %@", key);
//    }
    progressWindow.actionMessage.stringValue = @"Loading stack entries";
    gettimeofday(&time1, NULL);
//    [self loadStacksEntriesFromTagFile:stacksDocument];
    [self loadStacksEntriesFromTagFile:stacksDocument progressWindow:progressWindow];
    CHECK_STOP
    gettimeofday(&time2, NULL);
    NSLog(@"finished loading stacks entries time %ld", time2.tv_sec - time1.tv_sec);

    NSLog(@"loading snps onto datums");
    gettimeofday(&time1, NULL);
    progressWindow.actionMessage.stringValue = @"Loading datum snps";
    [self loadSnpsOntoDatum:stacksDocument];
    gettimeofday(&time2, NULL);
    NSLog(@"finished loading snps onto datum time %ld", time2.tv_sec - time1.tv_sec);

    NSLog(@"loading alleles onto datums");
    gettimeofday(&time1, NULL);
    progressWindow.actionMessage.stringValue = @"Loading datum alleles";
    [self loadAllelesOntoDatum:stacksDocument];
    gettimeofday(&time2, NULL);
    NSLog(@"finished loading alleles onto datum time %ld", time2.tv_sec - time1.tv_sec);


    NSError *error;
    progressWindow.actionMessage.stringValue = @"Final save";
    [stacksDocument.managedObjectContext save:&error];
//    NSLog(@"saved %d", success);
    [bar incrementBy:5];

    return stacksDocument;
}


- (void)loadSnpsOntoDatum:(StacksDocument *)document {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *path = document.path;

    // 1- get all files i the directory named *.tag.tsv"
    NSError *error;
    NSArray *files = [fileManager contentsOfDirectoryAtPath:path error:&error];
    NSLog(@"# of files for directory %ld", files.count);

    // 2 - for each file, read the .tags file
    for (NSString *filePath in files) {
        if ([filePath hasSuffix:@".snps.tsv"] && ![filePath hasPrefix:@"batch"]) {
            [self loadSnpFileForDatum:document fromFile:filePath];
        }
        else {
//            NSLog(@"not loading tag file %@", filePath);
        }
    }
}

- (void)loadAllelesOntoDatum:(StacksDocument *)document {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *path = document.path;

    // 1- get all files i the directory named *.tag.tsv"
    NSError *error;
    NSArray *files = [fileManager contentsOfDirectoryAtPath:path error:&error];
    NSLog(@"# of files for directory %ld", files.count);

    // 2 - for each file, read the .tags file
    for (NSString *filePath in files) {
        if ([filePath hasSuffix:@".alleles.tsv"] && ![filePath hasPrefix:@"batch"]) {
            [self loadAlleleFileForDatum:document fromFile:filePath];
        }
        else {
//            NSLog(@"not loading alleles file %@", filePath);
        }
    }
}

- (void)loadAlleleFileForDatum:(StacksDocument *)document fromFile:(NSString *)alleleFileName {
    NSLog(@"Loading allele file %@", alleleFileName);

    NSUInteger fileNameLength = alleleFileName.length;
    NSString *sampleName = [alleleFileName substringToIndex:fileNameLength - 12];
    NSLog(@"sampleName %@", sampleName);
    // sampleName . . . from lsat index of "/" . . . to just before ".tags.tsv"

    NSManagedObjectContext *moc = document.managedObjectContext;
    SampleMO *sampleMO = [[SampleRepository sharedInstance] getSampleForName:sampleName andContext:document.managedObjectContext andError:nil];

    struct timeval time1, time2;
    gettimeofday(&time1, NULL);

    NSError *error2 = nil;
    NSString *absoluteFileName = [document.path stringByAppendingFormat:@"/%@", alleleFileName];
    NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error2] componentsSeparatedByString:@"\n"];
    gettimeofday(&time2, NULL);
//    NSLog(@"load file data and split %ld", (time2.tv_sec - time1.tv_sec));

    if (error2) {
        NSLog(@"error loading file [%@]: %@", alleleFileName, error2);
    }

    NSString *line;
//    NSUInteger row = 1;
    gettimeofday(&time1, NULL);
    NSInteger locusId = -1;
    NSInteger newLocusId;
    DatumMO *datumMO = nil ;
    NSMutableDictionary *lookupDictionary = [sampleLookupDictionary objectForKey:sampleName];
//    char allele;
    int depth;
    float ratio;
//    LocusMO *locusMO = nil ;
    NSNumberFormatter *numberFormatter = [[NSNumberFormatter alloc] init];
    numberFormatter.numberStyle = NSNumberFormatterNoStyle;

    NSDictionary *datumLociMap = [[DatumRepository sharedInstance] getDatums:document.managedObjectContext forSample:sampleMO.sampleId];
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

        if (columns.count > 4) {

            // if the StackMO is found
//            sampleId = [[columns objectAtIndex:1] integerValue];
//            newLocusId = [[columns objectAtIndex:2] integerValue];

            NSString *internalIndex = (NSString *) [columns objectAtIndex:2];
            NSString *retrievedObject = (NSString *) [lookupDictionary objectForKey:internalIndex];
            newLocusId = [retrievedObject integerValue];


//            allele = [[columns objectAtIndex:3] charValue];
//            allele = [numberFormatter numberFromString:[columns objectAtIndex:3]];
            ratio = [[columns objectAtIndex:4] floatValue];
            depth = [[columns objectAtIndex:5] intValue];


            if (locusId != newLocusId) {
                locusId = newLocusId;
//                datumMO = [[DatumRepository sharedInstance] getDatum:moc locusId:locusId andSampleId:sampleMO.sampleId.integerValue];

                datumMO = [datumLociMap objectForKey:[NSString stringWithFormat:@"%ld", locusId]];
            }

            if (datumMO != nil) {
                // TODO: convert
//                [[AlleleRepository sharedInstance] insertDatumAllele:moc
//                                              ratio:[NSNumber numberWithFloat:ratio]
//                                              depth:[NSNumber numberWithInt:depth]
//                                             allele:[numberFormatter numberFromString:[columns objectAtIndex:3]]
//                                              datum:datumMO
//                ];
            }
        }
    }
    gettimeofday(&time2, NULL);


    // save old
    NSError *saveError;
    [moc save:&saveError];
//    NSLog(@"saved %d", success);
    if (saveError != nil ) {
        NSLog(@"error saving %@", saveError);
    }
}

- (void)loadSnpFileForDatum:(StacksDocument *)document fromFile:(NSString *)snpFileName {
    NSLog(@"Loading snps file %@", snpFileName);

    NSUInteger fileNameLength = snpFileName.length;
    NSString *sampleName = [snpFileName substringToIndex:fileNameLength - 9];
    NSLog(@"sampleName %@", sampleName);
    // sampleName . . . from lsat index of "/" . . . to just before ".tags.tsv"

    // create matches file name
    NSString *matchesFileName = [document.path stringByAppendingString:[sampleName stringByAppendingString:@".matches.tsv"]];
    NSLog(@"matches filename %@", matchesFileName);

//    NSMutableDictionary *lookupDictionary = [self loadMatchesDictionary:matchesFileName];
    NSMutableDictionary *lookupDictionary = [sampleLookupDictionary objectForKey:sampleName];
//    for(id key in [sampleLookupDictionary allKeys]){
//        NSLog(@"keys: %@",key);
//    }
    NSLog(@"size of lookupDictionary %ld", lookupDictionary.count);

    NSManagedObjectContext *moc = document.managedObjectContext;
    SampleMO *sampleMO = [[SampleRepository sharedInstance] getSampleForName:sampleName andContext:document.managedObjectContext andError:nil];

    struct timeval time1, time2;
    gettimeofday(&time1, NULL);

    NSError *error2 = nil;
    NSString *absoluteFileName = [document.path stringByAppendingFormat:@"/%@", snpFileName];
    NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error2] componentsSeparatedByString:@"\n"];
    gettimeofday(&time2, NULL);
//    NSLog(@"load file data and split %ld", (time2.tv_sec - time1.tv_sec));

    if (error2) {
        NSLog(@"error loading file [%@]: %@", snpFileName, error2);
    }

    NSString *line;
//    NSUInteger row = 1;
    gettimeofday(&time1, NULL);
    NSInteger locusId = -1;
//    NSInteger sampleId = -1;
    NSInteger newLocusId;
    NSInteger column;
    float lratio;
//    char rank1, rank2, rank3, rank4;
    DatumMO *datumMO = nil ;
//    LocusMO *locusMO = nil ;
    NSNumberFormatter *numberFormatter = [[NSNumberFormatter alloc] init];

    NSDictionary *datumLociMap = [[DatumRepository sharedInstance] getDatums:document.managedObjectContext forSample:sampleMO.sampleId];

    numberFormatter.numberStyle = NSNumberFormatterNoStyle;

    NSDictionary *snpLociMap = [self getLocusSnpsForDocument:document];
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

        if (columns.count > 6) {

            // if the StackMO is found
//            sampleId = [[columns objectAtIndex:1] integerValue];
//            newLocusId = [[columns objectAtIndex:2] integerValue];

            NSString *internalIndex = (NSString *) [columns objectAtIndex:2];
            NSString *retrievedObject = (NSString *) [lookupDictionary objectForKey:internalIndex];
//            NSLog(@"retrieved object: %@", retrievedObject);
            newLocusId = [retrievedObject integerValue];

            column = [[columns objectAtIndex:3] integerValue];
            lratio = [[columns objectAtIndex:4] floatValue];


            if (locusId != newLocusId) {
                locusId = newLocusId;
//                locusMO = [[LocusRepository sharedInstance] getLocus:moc forId:locusId];
                // search for the new locus
                // TODO: get from in-memory lookup?
//                datumMO = [[DatumRepository sharedInstance] getDatum:moc locusId:locusId andSampleId:sampleMO.sampleId.integerValue];
                datumMO = [datumLociMap objectForKey:[NSString stringWithFormat:@"%ld", locusId]];
            }

            if (datumMO != nil) {
                // TODO: convert
//                datumMO.snpData = [snpLociMap valueForKey:[NSString stringWithFormat:@"%ld",locusId]] ;

                // if we ahve data, then we want to add it


//                NSError *error;
//                NSDictionary *snpJson = [NSJSONSerialization JSONObjectWithData:datumMO.snpData options:kNilOptions error:&error];
//                NSMutableDictionary *snpLookup = [[NSMutableDictionary  alloc] init];
//                for (NSDictionary *snp in snpJson) {
////                    NSInteger startRange = [[snp valueForKey:@"column"] integerValue];
////                    [snpInts addObject:[snp valueForKey:@"column"]];
//                    [snpLookup setObject:snp forKey:[snp valueForKey:@"column"]];
////                    NSRange selectedRange = NSMakeRange(startRange, 1);
////                    [string addAttributes:attributes range:selectedRange];
//                }

                NSDictionary *snpDictionary = [NSDictionary dictionaryWithObjectsAndKeys:
                        [NSNumber numberWithInteger:column], @"column"
                        , [NSNumber numberWithFloat:lratio], @"lratio"
                        , [numberFormatter numberFromString:[columns objectAtIndex:5]], @"rank1"
                        , [numberFormatter numberFromString:[columns objectAtIndex:6]], @"rank2"
                        , [numberFormatter numberFromString:[columns objectAtIndex:7]], @"rank3"
                        , [numberFormatter numberFromString:[columns objectAtIndex:8]], @"rank4"
                        , nil ];
                
                NSMutableArray *snpArray;
                if(datumMO.snpData==nil){
                    snpArray = [[NSMutableArray alloc] init];
                }
                else{
                   snpArray = [NSMutableArray arrayWithArray:[NSJSONSerialization JSONObjectWithData:datumMO.snpData options:kNilOptions error:&error2]];
                }
                
                [snpArray addObject:snpDictionary];

                NSData *snpArrayData = [NSJSONSerialization dataWithJSONObject:snpArray options:NSJSONWritingPrettyPrinted error:&error2];
                datumMO.snpData = snpArrayData;

//                [[SnpRepository sharedInstance] insertDatumSnp:moc column:[NSNumber numberWithInteger:column]
//                                       lratio:[NSNumber numberWithFloat:lratio]
//                                        rank1:[numberFormatter numberFromString:[columns objectAtIndex:5]]
//                                        rank2:[numberFormatter numberFromString:[columns objectAtIndex:6]]
//                                        rank3:[numberFormatter numberFromString:[columns objectAtIndex:7]]
//                                        rank4:[numberFormatter numberFromString:[columns objectAtIndex:8]]
//                                        datum:datumMO
//                ];
//                [datumMO addSnpsObject:datumSnpMO];
//                NSLog(@"inserted snp at %@ for sample %@ and locus %@",datumSnpMO.column,datumMO.sample.name,datumMO.locus.locusId);
            }
        }
    }
    gettimeofday(&time2, NULL);


    // save old
    NSError *saveError;
    [moc save:&saveError];
    if (saveError != nil ) {
        NSLog(@"error saving %@", saveError);
    }

}

- (NSMutableDictionary *)loadMatchesDictionary:(NSString *)name {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    if (![fileManager fileExistsAtPath:name]) return nil;

//    NSString* contents = [NSString stringWithContentsOfFile:path encoding:NSUTF8StringEncoding]
    NSError *error2;
    NSArray *fileData = [[NSString stringWithContentsOfFile:name encoding:NSUTF8StringEncoding error:&error2] componentsSeparatedByString:@"\n"];
    NSMutableDictionary *lookupDictionary = [[NSMutableDictionary alloc] init];
    NSString *line;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];
        if (columns.count > 5) {
//        NSLog(@"column 2 values [%@] - [%@]",[columns objectAtIndex:2],[columns objectAtIndex:4]);
            NSString *internalId = [columns objectAtIndex:4];
            NSString *externalId = [columns objectAtIndex:2];
            [lookupDictionary setObject:externalId forKey:internalId];
        }
    }

    return lookupDictionary;
}

- (void)loadStacksEntriesFromTagFile:(StacksDocument *)document progressWindow:(ProgressController *)progressWindow {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *path = document.path;
    // TODO: if male . . .male.tags.tsv / female.tags.tsv
    // TODO: or *|!male|!female_N.tags.tsv

    // 1- get all files i the directory named *.tag.tsv"
    NSError *error;
//    NSArray *files = [fileManager contentsOfDirectoryAtPath:path error:&error];
    NSArray *files = [fileManager contentsOfDirectoryAtPath:path error:&error];
    NSMutableArray *realFiles = [[NSMutableArray alloc] init];

    for (NSString *filePath in files) {
        if ([filePath hasSuffix:@".tags.tsv"] && ![filePath hasPrefix:@"batch"]) {
            [realFiles addObject:filePath];
        }
    }

    NSLog(@"# of files for directory %ld", files.count);


    // 2 - for each file, read the .tags file
    int fileNumber = 0;
    NSUInteger numFiles = realFiles.count;
    double incrementAmount = 30.0 / numFiles;

    struct timeval time1, time2;
    gettimeofday(&time1, NULL);

    long saveAtLine = 500000;
    NSUInteger saveCounter = 1;

    NSMutableSet *savedDatums = [[NSMutableSet alloc] init];

    StackEntryRenderer *stackEntryView = [[StackEntryRenderer alloc] init];

    for (NSString *tagFileName in realFiles) {
        progressWindow.actionMessage.stringValue = [NSString stringWithFormat:@"Loading stack entry %i / %ld", fileNumber + 1, numFiles];
        if (stopProcess) return;
        if ([tagFileName hasSuffix:@".tags.tsv"] && ![tagFileName hasPrefix:@"batch"]) {
//            CHECK_STOP
//            [self loadTagFile:document fromFile:filePath];

            NSLog(@"Loading tag file %@", tagFileName);

            NSUInteger fileNameLength = tagFileName.length;
            NSString *sampleName = [tagFileName substringToIndex:fileNameLength - 9];
//    NSLog(@"sampleName %@", sampleName);
            // sampleName . . . from lsat index of "/" . . . to just before ".tags.tsv"

            NSManagedObjectContext *moc = document.managedObjectContext;
            SampleMO *sampleMO = [[SampleRepository sharedInstance] getSampleForName:sampleName andContext:document.managedObjectContext andError:nil];

            struct timeval time1, time2;
            gettimeofday(&time1, NULL);

            NSError *error2 = nil;
            NSString *absoluteFileName = [document.path stringByAppendingFormat:@"/%@", tagFileName];
            NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error2] componentsSeparatedByString:@"\n"];
//    NSLog(@"load file data and split %ld", (time2.tv_sec - time1.tv_sec));

            if (error2) {
                NSLog(@"error loading file [%@]: %@", tagFileName, error2);
            }

            NSString *line;
            NSUInteger row = 1;
            gettimeofday(&time1, NULL);
            NSInteger locusId = -1;
            NSInteger newLocusId;
            DatumMO *datumMO = nil ;

            NSDictionary *datumLociMap = [[DatumRepository sharedInstance] getDatums:document.managedObjectContext forSample:sampleMO.sampleId];
            NSDictionary *snpLociMap = [self getLocusSnpsForDocument:document];


            NSMutableDictionary *lookupDictionary = [sampleLookupDictionary objectForKey:sampleName];
//    for (id key in [sampleLookupDictionary allKeys]) {
//        NSLog(@"keys 3: %@", key);
//    }
//    NSLog(@"size of lookupDictionary %ld", lookupDictionary.count);


            for (line in fileData) {
                NSArray *columns = [line componentsSeparatedByString:@"\t"];

                if (columns.count > 8) {

                    // if the StackMO is found
//            newLocusId = [[columns objectAtIndex:2] integerValue];
                    NSString *internalIndex = (NSString *) [columns objectAtIndex:2];
                    NSString *retrievedObject = (NSString *) [lookupDictionary objectForKey:internalIndex];
//            NSLog(@"retrieved object: %@", retrievedObject);
                    newLocusId = [retrievedObject integerValue];
                    if (locusId != newLocusId) {

                        // cleanup old object if exists
//                datumMO.stackData = [NSString stringWithFormat:@"<p>Some stack data for sample '%@' and locus '%@'</p>",datumMO.sample.name,datumMO.locus.locusId];
                        if (false == [stackEntryView isEmpty]) {
                            stackEntryView.snpDatumData = datumMO.snpData ;
                            datumMO.stackData = [stackEntryView renderHtml];

//                    NSLog(@"saving for new locus");
//                    NSError *error;
//                    [moc save:&error];
//                    [moc refreshObject:datumMO mergeChanges:YES];

//                    datumMO = nil ;
//                    stackEntryView = nil ;
                        }


                        locusId = newLocusId;
                        // search for the new locus
                        // TODO: get from in-memory lookup?
//                datumMO = [[DatumRepository sharedInstance] getDatum:moc locusId:locusId andSampleId:sampleMO.sampleId.integerValue];
                        datumMO = [datumLociMap objectForKey:[NSString stringWithFormat:@"%ld", locusId]];

                        if (datumMO != nil) {
                            [savedDatums addObject:datumMO];
//                NSLog(@"%@ vs %@",sampleMO.sampleId,datumMO.sampleId );
//                        stackEntryView = [[StackEntryRenderer alloc] init];
                            [stackEntryView clear];

                            stackEntryView.locusId = locusId;

                            stackEntryView.snpLocusData = [snpLociMap objectForKey:[NSNumber numberWithInteger:locusId]];
//                            stackEntryView.sn= [snpLociMap objectForKey:[NSNumber numberWithInteger:locusId]];

//                        stackEntryView.sampleName = datumMO.sample.name;
                            stackEntryView.sampleName = sampleMO.name;

                            // TODO: map locus and datum snps
                        }
//                        else{
//                            NSLog(@"no datum found for locus: %ld and sample %@",locusId,sampleMO.name);
//                        }
                    }

                    if (datumMO != nil) {
                        NSString *relationship = [columns objectAtIndex:6];



                        if ([relationship isEqualToString:@"consensus"]) {
                            row = 1;
                            stackEntryView.consensus = [columns objectAtIndex:9];
//                            [stackEntryView.relationships addObject:@"consensus"];
//                    datumMO.consensus = [stackEntryRepository insertConsensusStackEntry:moc
//                                                                                  block:[columns objectAtIndex:7]
//                                                                             sequenceId:[columns objectAtIndex:8]
//                                                                               sequence:[columns objectAtIndex:9]
//                                                                                  datum:datumMO
//                    ];
//                    datumMO.reference = [stackEntryRepository insertReferenceStackEntry:moc
//                                                                               sequence:[columns objectAtIndex:9]
//                                                                                  datum:datumMO
//                    ];

                        }
                        else if ([relationship isEqualToString:@"model"]) {
                            stackEntryView.model = [columns objectAtIndex:9];
//                            [stackEntryView.relationships addObject:@"model"];
//                    datumMO.model = [stackEntryRepository insertModelStackEntry:moc
//                                                                          block:[columns objectAtIndex:7]
//                                                                     sequenceId:[columns objectAtIndex:8]
//                                                                       sequence:[columns objectAtIndex:9]
//                                                                          datum:datumMO
//                    ];
                        }
                        else {
                            [stackEntryView.sequenceIds addObject:[columns objectAtIndex:8]];
                            [stackEntryView.sequences addObject:[columns objectAtIndex:9]];
                            [stackEntryView.relationships addObject:[columns objectAtIndex:6]];
                            [stackEntryView.blocks addObject:[columns objectAtIndex:7]];
                            [stackEntryView.entryIds addObject:[NSNumber numberWithInteger:row]];
//                    StackEntryMO *stackEntryMO = [stackEntryRepository insertStackEntry:moc
//                                                                                entryId:[NSNumber numberWithInteger:row]
//                                                                           relationship:[columns objectAtIndex:6]
//                                                                                  block:[columns objectAtIndex:7]
//                                                                             sequenceId:[columns objectAtIndex:8]
//                                                                               sequence:[columns objectAtIndex:9]
//                                                                              consensus:datumMO.consensus.sequence
//                                                                                  datum:datumMO
//                    ];
//                    [datumMO addStackEntriesObject:stackEntryMO];
                            ++row;
                        }

//                datumMO.stackData = [NSString stringWithFormat:@"<p>Some stack data for sample '%@' and locus '%@'</p>",datumMO.sample.name,datumMO.locus.locusId];


                        if (saveCounter % saveAtLine == 0) {
                            NSLog(@"SAVING");
                            NSError *saveError;
                            [moc save:&saveError];
                            if (saveError != nil ) {
                                NSLog(@"error saving %@", saveError);
                            }
                            for (DatumMO *datumMO in savedDatums) {
                                [moc refreshObject:datumMO mergeChanges:YES];
                            }
                            [savedDatums removeAllObjects];
                            NSLog(@"saved datums count: %ld", savedDatums.count);
                        }


                        ++saveCounter;
                    }


                }

            }



//            if (saveCounter % saveAtLine == 0) {
//                NSLog(@"SAVING");
//                NSError *saveError;
//                [moc save:&saveError];
//                if (saveError != nil ) {
//                    NSLog(@"error saving %@", saveError);
//                }
//            }


            gettimeofday(&time2, NULL);
//    NSLog(@"parse entries lines %ld produce %ld - %ld", fileData.count, datumMO.stackEntries.count, (time2.tv_sec - time1.tv_sec));
            NSLog(@"parse entries lines %ld time: %ld", fileData.count, (time2.tv_sec - time1.tv_sec));


        }
        else {
            NSLog(@"not loading tag file %@", tagFileName);
        }
        [progressWindow.loadProgress incrementBy:incrementAmount];
        ++fileNumber;
    }
    gettimeofday(&time2, NULL);

    NSLog(@"STACK ENTRY LOAD: time for %ld  and loci %ld file: %ld s", numFiles, document.loci.count, (time2.tv_sec - time1.tv_sec));

    // save old
    NSError *saveError;
    [document.managedObjectContext save:&saveError];
    NSLog(@"regular save");
    if (saveError != nil ) {
        NSLog(@"error saving %@", saveError);
    }

}

- (NSDictionary *)getLocusSnpsForDocument:(StacksDocument *)document {
    if(locusSnpMap==nil){
        locusSnpMap = [[NSMutableDictionary alloc] init];
        for(LocusMO *locusMO in document.loci){
            [locusSnpMap setObject:locusMO.snpData forKey:locusMO.locusId];
        }
    }
    return locusSnpMap;
}

//- (void)loadTagFile:(StacksDocument *)document fromFile:(NSString *)tagFileName {
//    NSLog(@"Loading tag file %@", tagFileName);
//
//    NSUInteger fileNameLength = tagFileName.length;
//    NSString *sampleName = [tagFileName substringToIndex:fileNameLength - 9];
////    NSLog(@"sampleName %@", sampleName);
//    // sampleName . . . from lsat index of "/" . . . to just before ".tags.tsv"
//
//    NSManagedObjectContext *moc = document.managedObjectContext;
//    SampleMO *sampleMO = [[SampleRepository sharedInstance] getSampleForName:sampleName andContext:document.managedObjectContext andError:nil];
//
//    struct timeval time1, time2;
//    gettimeofday(&time1, NULL);
//
//    NSError *error2 = nil;
//    NSString *absoluteFileName = [document.path stringByAppendingFormat:@"/%@", tagFileName];
//    NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error2] componentsSeparatedByString:@"\n"];
////    NSLog(@"load file data and split %ld", (time2.tv_sec - time1.tv_sec));
//
//    if (error2) {
//        NSLog(@"error loading file [%@]: %@", tagFileName, error2);
//    }
//
//    NSString *line;
//    NSUInteger row = 1;
//    gettimeofday(&time1, NULL);
//    NSInteger locusId = -1;
//    NSInteger newLocusId;
//    DatumMO *datumMO = nil ;
//
//    NSDictionary *datumLociMap = [[DatumRepository sharedInstance] getDatums:document.managedObjectContext forSample:sampleMO.sampleId];
//
//    StackEntryRenderer *stackEntryView = nil ;
//
//    NSMutableDictionary *lookupDictionary = [sampleLookupDictionary objectForKey:sampleName];
////    for (id key in [sampleLookupDictionary allKeys]) {
////        NSLog(@"keys 3: %@", key);
////    }
////    NSLog(@"size of lookupDictionary %ld", lookupDictionary.count);
//
//    NSUInteger saveAtLine = 50000;
//    NSUInteger saveCounter = 1;
//
//    for (line in fileData) {
//        NSArray *columns = [line componentsSeparatedByString:@"\t"];
//
//        if (columns.count > 8) {
//
//            // if the StackMO is found
////            newLocusId = [[columns objectAtIndex:2] integerValue];
//            NSString *internalIndex = (NSString *) [columns objectAtIndex:2];
//            NSString *retrievedObject = (NSString *) [lookupDictionary objectForKey:internalIndex];
////            NSLog(@"retrieved object: %@", retrievedObject);
//            newLocusId = [retrievedObject integerValue];
//            if (locusId != newLocusId) {
//
//                // cleanup old object if exists
////                datumMO.stackData = [NSString stringWithFormat:@"<p>Some stack data for sample '%@' and locus '%@'</p>",datumMO.sample.name,datumMO.locus.locusId];
//                if (stackEntryView != nil) {
//                    datumMO.stackData = [stackEntryView renderHtml];
//
////                    NSLog(@"saving for new locus");
////                    NSError *error;
////                    [moc save:&error];
////                    [moc refreshObject:datumMO mergeChanges:YES];
//
////                    datumMO = nil ;
////                    stackEntryView = nil ;
//                }
//
//
//                locusId = newLocusId;
//                // search for the new locus
//                // TODO: get from in-memory lookup?
////                datumMO = [[DatumRepository sharedInstance] getDatum:moc locusId:locusId andSampleId:sampleMO.sampleId.integerValue];
//                datumMO = [datumLociMap objectForKey:[NSString stringWithFormat:@"%ld", locusId]];
////                NSLog(@"%@ vs %@",sampleMO.sampleId,datumMO.sampleId );
//                stackEntryView = [[StackEntryRenderer alloc] init];
//                stackEntryView.locusId = locusId;
//                stackEntryView.sampleName = datumMO.sample.name;
//
//                // TODO: map locus and datum snps
//            }
//
//            if (datumMO != nil) {
//                NSString *relationship = [columns objectAtIndex:6];
//
//                if ([relationship isEqualToString:@"consensus"]) {
//                    row = 1;
//                    stackEntryView.consensus = [columns objectAtIndex:9];
////                    datumMO.consensus = [stackEntryRepository insertConsensusStackEntry:moc
////                                                                                  block:[columns objectAtIndex:7]
////                                                                             sequenceId:[columns objectAtIndex:8]
////                                                                               sequence:[columns objectAtIndex:9]
////                                                                                  datum:datumMO
////                    ];
////                    datumMO.reference = [stackEntryRepository insertReferenceStackEntry:moc
////                                                                               sequence:[columns objectAtIndex:9]
////                                                                                  datum:datumMO
////                    ];
//
//                }
//                else if ([relationship isEqualToString:@"model"]) {
//                    stackEntryView.model = [columns objectAtIndex:9];
////                    datumMO.model = [stackEntryRepository insertModelStackEntry:moc
////                                                                          block:[columns objectAtIndex:7]
////                                                                     sequenceId:[columns objectAtIndex:8]
////                                                                       sequence:[columns objectAtIndex:9]
////                                                                          datum:datumMO
////                    ];
//                }
//                else {
//                    [stackEntryView.sequenceIds addObject:[columns objectAtIndex:8]];
//                    [stackEntryView.sequences addObject:[columns objectAtIndex:9]];
////                    StackEntryMO *stackEntryMO = [stackEntryRepository insertStackEntry:moc
////                                                                                entryId:[NSNumber numberWithInteger:row]
////                                                                           relationship:[columns objectAtIndex:6]
////                                                                                  block:[columns objectAtIndex:7]
////                                                                             sequenceId:[columns objectAtIndex:8]
////                                                                               sequence:[columns objectAtIndex:9]
////                                                                              consensus:datumMO.consensus.sequence
////                                                                                  datum:datumMO
////                    ];
////                    [datumMO addStackEntriesObject:stackEntryMO];
//                    ++row;
//                }
//
////                datumMO.stackData = [NSString stringWithFormat:@"<p>Some stack data for sample '%@' and locus '%@'</p>",datumMO.sample.name,datumMO.locus.locusId];
//
//
//                if (saveCounter % saveAtLine == 0) {
//                    NSLog(@"SAVING");
//                    NSError *saveError;
//                    [moc save:&saveError];
//                    if (saveError != nil ) {
//                        NSLog(@"error saving %@", saveError);
//                    }
//                }
//
//
//                ++saveCounter;
//            }
//
//
//        }
//
//    }
//
//    if (saveCounter % saveAtLine == 0) {
//        NSLog(@"SAVING");
//        NSError *saveError;
//        [moc save:&saveError];
//        if (saveError != nil ) {
//            NSLog(@"error saving %@", saveError);
//        }
//    }
//
//
//    gettimeofday(&time2, NULL);
////    NSLog(@"parse entries lines %ld produce %ld - %ld", fileData.count, datumMO.stackEntries.count, (time2.tv_sec - time1.tv_sec));
//    NSLog(@"parse entries lines %ld time: %ld", fileData.count, (time2.tv_sec - time1.tv_sec));
//
//    datumMO = nil ;
//    stackEntryView = nil ;
//
////    // save old
////    NSError *saveError;
////    [moc save:&saveError];
////    NSLog(@"regular save") ;
////    if (saveError != nil ) {
////        NSLog(@"error saving %@", saveError);
////    }
////    [moc reset];
//}

- (void)addSamplesToDocument:(StacksDocument *)document forSampleIds:(vector<int>)sampleIds andSamples:(map<int, string>)samples {

    if (document.populationLookup == nil || document.populationLookup.count == 0) {
        PopulationMO *populationMO = [[PopulationRepository sharedInstance] insertPopulation:document.managedObjectContext id:[NSNumber numberWithInt:1] name:@"All"];
        document.populations = [NSSet setWithObjects:populationMO, nil];

        // set each sample to populationMO
        for (int i = 0; i < sampleIds.size(); i++) {
            SampleMO *sampleMO = [[SampleRepository sharedInstance] insertSample:document.managedObjectContext id:[NSNumber numberWithInt:sampleIds[i]] name:[NSString stringWithUTF8String:samples[sampleIds[i]].c_str()]];

            [populationMO addSamplesObject:sampleMO];
        }
    }
    else {
        NSNumberFormatter *f = [[NSNumberFormatter alloc] init];
        // else
        for (int i = 0; i < sampleIds.size(); i++) {
            SampleMO *sampleMO = [[SampleRepository sharedInstance] insertSample:document.managedObjectContext
                                                                              id:[NSNumber numberWithInt:sampleIds[i]]
                                                                            name:[NSString stringWithUTF8String:samples[sampleIds[i]].c_str()]];

            NSString *populationId = [document.populationLookup objectForKey:sampleMO.name];

            if (populationId != nil) {
                // lets get the population . . can use lookup, but this is usually pretty small
                NSLog(@"tyring to populate popid %@", populationId);
                for (PopulationMO *populationMO in document.populations) {
                    NSNumber *endNumber = [f numberFromString:populationId];
                    NSLog(@"comparing to %@ vs %@", populationMO.populationId, endNumber);
                    if ([populationMO.populationId isEqualToNumber:endNumber]) {
                        NSLog(@"FOuND population ID %@", populationMO.populationId);
                        [populationMO addSamplesObject:sampleMO];
                    }
                }
            }
        }
    }

}

- (void)addPopulationsToDocument:(StacksDocument *)document forPath:(NSString *)path {
    // sample . . . index / name
    NSMutableDictionary *populationLookup = [self loadPopulation:path];

    document.populationLookup = populationLookup;

    NSLog(@"population lookup %ld", document.populationLookup.count);


    // scan and create the populations . . .
    NSNumberFormatter *f = [[NSNumberFormatter alloc] init];
//    [f setNumberStyle:NSNumberFormatterDecimalStyle];

    // creates a unique set
    NSMutableSet *populationIdsSet = [[NSMutableSet alloc] init];
    for (NSString *populationId in populationLookup.allValues) {
        [populationIdsSet addObject:populationId];
    }

    for (NSString *popId in populationIdsSet) {
        NSLog(@"adding population %@", popId);
        NSNumber *myNumber = [f numberFromString:popId];
        [[PopulationRepository sharedInstance] insertPopulation:document.managedObjectContext id:myNumber name:popId];
    }

    NSArray *populationArray = [[PopulationRepository sharedInstance] getAllPopulations:document.managedObjectContext];

    document.populations = [NSSet setWithArray:populationArray];

}

- (void)readPopulations:(StacksDocument *)document {
    NSMutableArray *populations = [[NSMutableArray alloc] init];
    NSManagedObjectContext *moc = document.managedObjectContext;
    [moc setUndoManager:nil];

    NSString *path = document.path;
    NSLog(@"reading population for file %@", path);

    NSString *popmapFile = [path stringByAppendingString:@"popmap"];
    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL exists = [fileManager fileExistsAtPath:popmapFile];

    if (exists) {
        NSArray *fileData = [[NSString stringWithContentsOfFile:popmapFile encoding:NSUTF8StringEncoding error:nil] componentsSeparatedByString:@"\n"];
        NSString *line;
        for (line in fileData) {
            NSArray *columns = [line componentsSeparatedByString:@"\t"];
            if (columns.count == 2) {

//                [populationLookup setObject:[columns objectAtIndex:1] forKey:[columns objectAtIndex:0]];
//                NSString *sampleName = [columns objectAtIndex:0]; // the sample name . . . male, female, progeny, etc.

                NSString *populationName = [columns objectAtIndex:1]; // initially an integer
                PopulationMO *newPopulationMO = [[PopulationRepository sharedInstance] getPopulation:moc name:populationName];

                if (newPopulationMO != nil) {
                    [populations addObject:newPopulationMO];
                }
            }
            else {
                NSLog(@"something wrong with the column count %@", line);
            }
        }
    }
    else {
        NSLog(@"does not exist at %@", popmapFile);
    }

    NSLog(@"population size %ld", populations.count);
//    [document.populations setWithArray:populations ];
    document.populations = [[NSSet alloc] initWithArray:populations];
//    return populationLookup;

}

- (StacksDocument *)createStacksDocumentForPath:(NSString *)path {
    NSError *stacksDocumentCreateError;
    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
    stacksDocument.path = path;
    stacksDocument.name = path.lastPathComponent;
//    NSManagedObjectContext *context = [[NSManagedObjectContext alloc] initWithConcurrencyType:NSPrivateQueueConcurrencyType];

//    NSManagedObjectContext *moc = [stacksDocument getContextForPath:path];
    NSManagedObjectContext *moc = [self getContextForPath:path andName:path.lastPathComponent andDocument:stacksDocument];
    [moc setUndoManager:nil];
    stacksDocument.managedObjectContext = moc;
    if (stacksDocumentCreateError) {
        NSLog(@"error creating stacks document %@", stacksDocumentCreateError);
        return nil;
    }
    return stacksDocument;
}

- (void)setParentCounts:(NSManagedObjectContext *)context forFile:(NSString *)file {
    NSArray *fileData = [[NSString stringWithContentsOfFile:file encoding:NSUTF8StringEncoding error:nil] componentsSeparatedByString:@"\n"];
    NSString *line;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];
// should be column 8
        if (columns.count > 9) {
            NSInteger locusId = [[NSString stringWithFormat:@"%@", columns[2]] integerValue];
//            NSNumber *lookupKey = [NSNumber numberWithInteger:[[NSString stringWithFormat:@"%d", it->first] integerValue]];
            NSArray *parents = [columns[8] componentsSeparatedByString:@","];

            NSInteger parentCount = countParents(parents);
//            NSLog(@"parent count for %ld is %ld",locusId,parentCount);

            LocusMO *locusMO = [[LocusRepository sharedInstance] getLocus:context forId:locusId];
            locusMO.parentCount = [NSNumber numberWithInteger:parentCount];
        }
    }
}

- (void)dealloc {
    [sampleLookupDictionary removeAllObjects];
    persistentStoreCoordinator = nil ;
}

//- (StacksDocument *)getStacksDocumentForPath:(NSString *)path {
//    NSError *stacksDocumentCreateError;
//    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
//    NSManagedObjectContext *moc = [stacksDocument getContextForPath:path];
////    stacksDocument.managedObjectContext = moc ;
////    stacksDocument.path = path;
//    if (stacksDocumentCreateError) {
//        NSLog(@"error creating stacks document %@", stacksDocumentCreateError);
//        return nil;
//    }
//    return stacksDocument;
//}
@end

NSUInteger countParents(NSArray *parents) {
    NSUInteger parentCount = 0;

    NSMutableSet *parentIds = [[NSMutableSet alloc] initWithCapacity:parents.count];
    for (NSString *parentString in parents) {
        [parentIds addObject:[parentString componentsSeparatedByString:@"_"][0]];
    }
    parentCount = parentIds.count;

    return parentCount;
}



//void setParentCounts(NSString *file) {
//    NSArray *fileData = [[NSString stringWithContentsOfFile:file encoding:NSUTF8StringEncoding error:nil] componentsSeparatedByString:@"\n"];
//    NSString *line;
//    for (line in fileData) {
//        NSArray *columns = [line componentsSeparatedByString:@"\t"];
//        // should be column 8
//        if (columns.count > 9) {
//            NSNumber *locusId = [NSNumber numberWithInteger:[[NSString stringWithFormat:@"%@", columns[2]] integerValue]];
////            NSNumber *lookupKey = [NSNumber numberWithInteger:[[NSString stringWithFormat:@"%d", it->first] integerValue]];
//            NSArray *parents = [columns[8] componentsSeparatedByString:@","];
//
//            NSInteger parentCount = countParents(parents);
//            NSLog(@"parent count for %@ is %ld",locusId,parentCount);
//
//            ((LocusMO *) [dictionary objectForKey:locusId]).parentCount = [NSNumber numberWithUnsignedInt:parentCount];
//        }
//    }
//}

//void setParentCounts(NSMutableDictionary *dictionary, NSString *file) {
//    NSArray *fileData = [[NSString stringWithContentsOfFile:file encoding:NSUTF8StringEncoding error:nil] componentsSeparatedByString:@"\n"];
//    NSString *line;
//    for (line in fileData) {
//        NSArray *columns = [line componentsSeparatedByString:@"\t"];
//        // should be column 8
//        if (columns.count > 9) {
//            NSNumber *locusId = [NSNumber numberWithInteger:[[NSString stringWithFormat:@"%@", columns[2]] integerValue]];
////            NSNumber *lookupKey = [NSNumber numberWithInteger:[[NSString stringWithFormat:@"%d", it->first] integerValue]];
//            NSArray *parents = [columns[8] componentsSeparatedByString:@","];
//
//            NSInteger parentCount = countParents(parents);
//            NSLog(@"parent count for %@ is %ld",locusId,parentCount);
//            ((LocusMO *) [dictionary objectForKey:locusId]).parentCount = [NSNumber numberWithUnsignedInt:parentCount];
//        }
//    }
//}


NSString *calculateType(NSString *file) {
    NSArray *fileData = [[NSString stringWithContentsOfFile:file encoding:NSUTF8StringEncoding error:nil] componentsSeparatedByString:@"\n"];
    NSString *line;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];
        // should be column 8
//        NSLog(@"count is %ld",columns.count);
        if (columns.count > 9) {
            NSArray *parents = [columns[8] componentsSeparatedByString:@","];
            NSInteger parentCount = countParents(parents);
//            NSLog(@"parents %@ count %ld",columns[8],parentCount);
            if (parentCount > 2) {
                return @"Population";
            }
        }
    }
    return @"GeneticMap";
}


