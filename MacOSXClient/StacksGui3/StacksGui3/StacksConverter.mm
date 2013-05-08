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


#import "StackEntryMO.h"
#import "SampleRepository.h"
// #import "GenotypeView.h"
//#import "LocusView.h"
#import "StacksDocument.h"
#import "StacksConverter.h"
#import "PopulationRepository.h"
#import "DepthRepository.h"
#import "DatumRepository.h"
#import "HaplotypeRepository.h"
#import "LocusRepository.h"
//#import "StacksView.h"
// #import "DataStubber.h"


// #include "LociLoader.hpp"
//#import "GenotypeEntry.h"
//#import "SnpView.h"
#import "SnpMO.h"
#import "DatumMO.h"
#import "HaplotypeMO.h"
#import "DepthMO.h"
#import "PopulationMO.h"
#import "SampleMO.h"
#import "StackMO.h"
#import "SnpRepository.h"
#import "StackEntryRepository.h"
#import "StackRepository.h"
//#import "StackEntry.h"


#include <sys/time.h>

@implementation StacksConverter {

}

@synthesize datumRepository;
@synthesize depthRepository;
@synthesize haplotypeRepository;
@synthesize locusRepository;
@synthesize populationRepository;
@synthesize sampleRepository;
@synthesize snpRepository;
@synthesize stackEntryRepository;
@synthesize stackRepository;


- (id)init {
    self = [super init];
    if (self) {
        // nothing write now
        datumRepository = [[DatumRepository alloc] init];
        depthRepository = [[DepthRepository alloc] init];
        haplotypeRepository = [[HaplotypeRepository alloc] init];
        locusRepository = [[LocusRepository alloc] init];
        populationRepository = [[PopulationRepository alloc] init];
        sampleRepository = [[SampleRepository alloc] init];
        snpRepository = [[SnpRepository alloc] init];
        stackEntryRepository = [[StackEntryRepository alloc] init];
        stackRepository = [[StackRepository alloc] init];
    }
    return self;
}

//- (StacksView *)loadStacksView:(NSString *)filename atPath:(NSString *)path forTag:(NSInteger)tag locus:(LocusView *)locus {
//    NSFileManager *fileManager = [NSFileManager defaultManager];
//    // TODO: if male . . .male.tags.tsv / female.tags.tsv
//    // TODO: or *|!male|!female_N.tags.tsv
//
//    NSString *absoluteFileName = [NSString stringWithFormat:@"%@/%@.tags.tsv", path, filename];
//
//    BOOL existsAtPath = [fileManager fileExistsAtPath:absoluteFileName];
//    NSLog(@"absolute file name %@ exists %d", absoluteFileName, existsAtPath);
//    if (!existsAtPath) {
//        return nil;
//    }
//    StacksView *stacksView = [[StacksView alloc] init];
//
//
//    struct timeval time1, time2;
//    gettimeofday(&time1, NULL);
//
//    NSError *error = nil;
//    NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error] componentsSeparatedByString:@"\n"];
//    gettimeofday(&time2, NULL);
//    NSLog(@"load file data and split %ld", (time2.tv_sec - time1.tv_sec));
//
//    if (error) {
//        NSLog(@"error loading file [%@]: %@", absoluteFileName, error);
//    }
//
//    NSMutableArray *stackEntries = [[NSMutableArray alloc] init];
//
//    NSString *line;
//    NSUInteger row = 1;
//    gettimeofday(&time1, NULL);
//    for (line in fileData) {
//        NSArray *columns = [line componentsSeparatedByString:@"\t"];
//
//        if (columns.count > 8 && [[columns objectAtIndex:2] integerValue] == tag) {
//            StackEntryMO *stackEntry = [[StackEntryMO alloc] init];
////            stackEntry.entryId = row;
//            stackEntry.entryId = [NSNumber numberWithInteger:row];
//            stackEntry.relationship = [columns objectAtIndex:6];
//            stackEntry.block = [columns objectAtIndex:7];
//            stackEntry.sequenceId = [columns objectAtIndex:8];
//            stackEntry.sequence = [columns objectAtIndex:9];
//
//            if ([stackEntry.relationship isEqualToString:@"consensus"]) {
//                stacksView.consensus = stackEntry;
//            }
//            else if ([stackEntry.relationship isEqualToString:@"model"]) {
//                stacksView.model = stackEntry;
//            }
//            else {
//                ++row;
//                [stackEntries addObject:stackEntry];
//            }
//        }
//    }
//    gettimeofday(&time2, NULL);
//    NSLog(@"parse entries lines %ld produce %ld - %ld", fileData.count, stackEntries.count, (time2.tv_sec - time1.tv_sec));
//
//    stacksView.stackEntries = stackEntries;
//    StackEntryMO *referenceStack = [[StackEntryMO alloc] init];
//    referenceStack.entryId = nil ;
//    stacksView.reference = referenceStack;
//
//    // random snps
//    NSMutableArray *snps = [[NSMutableArray alloc] init];
//
//    for (SnpView *snp in locus.snps) {
//        [snps addObject:[NSNumber numberWithInt:snp.column]];
//    }
//
//    stacksView.snps = snps;
//
//    return stacksView;
//}

- (StacksDocument *)loadLociAndGenotypes:(NSString *)path {

    StacksDocument *stacksDocument = [self createStacksDocumentForPath:path];
    if (stacksDocument == nil) {
        return nil;
    }

    return [self loadDocument:stacksDocument];

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
- (void)checkFile:(NSString *)examplePath {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];

    if (!existsAtPath) {
        NSLog(@"files do not exist %@", examplePath);
        exit(1);
    }
}

- (StacksDocument *)loadDocument:(StacksDocument *)stacksDocument {
    NSString *path = stacksDocument.path;
    [self checkFile:path];
    map<int, CSLocus *> catalog;
    NSString *catalogFile = [path stringByAppendingString:@"batch_1.catalog"];

    struct timeval time1, time2;
    gettimeofday(&time1, NULL);


    load_loci([catalogFile UTF8String], catalog, false);
    gettimeofday(&time2, NULL);
    NSLog(@"load_loci %ld", (time2.tv_sec - time1.tv_sec));



    /**
    * START: loading datums + extra info
*/
    vector<vector<CatMatch *> > catalog_matches;
    map<int, string> samples;
    vector<int> sample_ids;

    vector<pair<int, string>> files = [self buildFileList:path];
    NSLog(@"number of files %ld", files.size());



    // loci loaded . . . now loading datum
    NSLog(@"model size %d", (int) catalog.size());

    gettimeofday(&time1, NULL);
    for (uint i = 0; i < files.size(); i++) {
        vector<CatMatch *> m;
        load_catalog_matches([[path stringByAppendingString:@"/"] UTF8String] + files[i].second, m);

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
    }
    gettimeofday(&time2, NULL);
    NSLog(@"catalog matches %ld", (time2.tv_sec - time1.tv_sec));


    NSLog(@"input populations %ld", stacksDocument.populations.count);
    [self addPopulationsToDocument:stacksDocument forPath:path];
    NSLog(@"output populations %ld", stacksDocument.populations.count);
    [self addSamplesToDocument:stacksDocument forSampleIds:sample_ids andSamples:samples];
    NSLog(@"output populations after samples %ld", stacksDocument.populations.count);
    for (PopulationMO *populationMo in stacksDocument.populations) {
        NSLog(@"samples %ld per population %@", populationMo.samples.count, populationMo.name);
//        for(SampleMO* sampleMO in populationMo.samples){
//            NSLog(@"sample %@ ",sampleMO.name);
//        }
    }

    cerr << "Populating observed haplotypes for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";

    gettimeofday(&time1, NULL);
    PopMap<CSLocus> *pmap = new PopMap<CSLocus>(sample_ids.size(), catalog.size());
    pmap->populate(sample_ids, catalog, catalog_matches);
    gettimeofday(&time2, NULL);
    NSLog(@"population pmap %ld", (time2.tv_sec - time1.tv_sec));


    NSMutableSet *loci = [[NSMutableSet alloc] init];

    gettimeofday(&time1, NULL);
    map<int, CSLocus *>::iterator catalogIterator = catalog.begin();
    NSManagedObjectContext *moc = stacksDocument.managedObjectContext;
    while (catalogIterator != catalog.end()) {
        const char *read = (*catalogIterator).second->con;
        LocusMO *locusMO = [locusRepository insertNewLocus:moc withId:[NSNumber numberWithInt:(*catalogIterator).first]
                                              andConsensus:[[NSString alloc] initWithCString:read encoding:NSUTF8StringEncoding] andMarker:[NSString stringWithUTF8String:catalogIterator->second->marker.c_str()]
        ];
        vector<SNP *> snps = catalogIterator->second->snps;
        vector<SNP *>::iterator snpsIterator = snps.begin();

        for (; snpsIterator != snps.end(); ++snpsIterator) {
            SNP *snp = (*snpsIterator);

            SnpMO *snpMO = [snpRepository insertSnp:moc
                                             column:[NSNumber numberWithInt:snp->col]
                                             lratio:[NSNumber numberWithFloat:snp->lratio]
                                              rank1:[NSNumber numberWithInt:snp->rank_1]
                                              rank2:[NSNumber numberWithInt:snp->rank_2]
                                              rank3:[NSNumber numberWithInt:snp->rank_3]
                                              rank4:[NSNumber numberWithInt:snp->rank_4]
            ];
            [locusMO addSnpsObject:snpMO];
        }

        [loci addObject:locusMO];
        ++catalogIterator;
    }
    gettimeofday(&time2, NULL);
    NSLog(@"populating snps %ld", (time2.tv_sec - time1.tv_sec));

    map<int, CSLocus *>::iterator it;
    Datum *datum;
    CSLocus *loc;

    // for each sample process the catalog
    NSLog(@"samples %ld X catalog %ld = %ld ", sample_ids.size(), catalog.size(), sample_ids.size() * catalog.size());

    long totalCatalogTime = 0;
    //go through all samples
    for (uint i = 0; i < sample_ids.size(); i++) {
        int sampleId = sample_ids[i];
        string sampleString = samples[sampleId];


        gettimeofday(&time1, NULL);
        // go through all loci
        for (it = catalog.begin(); it != catalog.end(); it++) {
            loc = it->second;
            datum = pmap->datum(loc->id, sample_ids[i]);

            LocusMO *locusMO = nil ;
            NSArray *locusArray = [loci allObjects];
            for (LocusMO *aLocus in locusArray) {
                NSNumber *lookupKey = [NSNumber numberWithInteger:[[NSString stringWithFormat:@"%ld", it->first] integerValue]];
                if ([lookupKey isEqualToNumber:aLocus.locusId]) {
                    locusMO = aLocus;
                }
            }

            if (datum != NULL && locusMO != nil) {
                NSString *key = [NSString stringWithUTF8String:sampleString.c_str()];


                NSError *error;
                SampleMO *sampleMO = [sampleRepository getSampleForName:key andContext:moc andError:nil];

                vector<char *> obshape = datum->obshap;
                vector<int> depths = datum->depth;
                int numLetters = obshape.size();
                DatumMO* newDatumMO = [datumRepository insertDatum:moc name:key sampleId:[NSNumber numberWithInt:sample_ids[i]] sample:sampleMO] ;

                // get catalogs for matches
                newDatumMO.tagId = [NSNumber numberWithInt:datum->id];
                locusMO.length = [NSNumber numberWithInt:loc->depth];

                if (depths.size() == numLetters) {
                    for (int j = 0; j < numLetters; j++) {
                        HaplotypeMO *haplotypeMO = [haplotypeRepository insertHaplotype:moc haplotype:[NSString stringWithUTF8String:obshape[j]]] ;
                        [newDatumMO addHaplotypesObject:haplotypeMO];

                        DepthMO *depthMO = [depthRepository insertDepth:moc depth:[NSNumber numberWithInt:depths[j]]];
                        [newDatumMO addDepthsObject:depthMO];
                    }
                    [locusMO addDatumsObject:newDatumMO];
                }
                else {
                    NSLog(@"mismatchon %@", [NSString stringWithUTF8String:sampleString.c_str()]);
                }

                vector<SNP *> snps = datum->snps;
                vector<SNP *>::iterator snpsIterator = snps.begin();

                if (snps.size() > 0) {
                    NSLog(@"has snps %ld", snps.size());
                }

//        NSMutableSet *snpsSet = [[NSMutableSet alloc] init];
                for (; snpsIterator != snps.end(); ++snpsIterator) {
                    SnpMO *snpMO = [NSEntityDescription insertNewObjectForEntityForName:@"Snp" inManagedObjectContext:stacksDocument.managedObjectContext];
                    SNP *snp = (*snpsIterator);
                    snpMO.column = [NSNumber numberWithInt:snp->col];
                    snpMO.lratio = [NSNumber numberWithFloat:snp->lratio];
                    snpMO.rank1 = [NSNumber numberWithChar:snp->rank_1];
                    snpMO.rank2 = [NSNumber numberWithChar:snp->rank_2];
                    snpMO.rank3 = [NSNumber numberWithChar:snp->rank_3];
                    snpMO.rank4 = [NSNumber numberWithChar:snp->rank_4];
                    [newDatumMO addSnpsObject:snpMO];
                }

//                }
//                else {
//                    NSLog(@"datum %@ FOUND for key %@ and locus %@", datumMO.name, key, locusMO.locusId);
//                }
//                datumMO.sample = sampleMO;
//                NSp
//                NSLog(@"sampleMO %@",sampleMO) ;
//                [sampleMO addDatumsObject:datumMO];
            }
        }
        gettimeofday(&time2, NULL);
        NSLog(@"iterating sample %d - time %ld", sample_ids[i], (time2.tv_sec - time1.tv_sec));
        totalCatalogTime += time2.tv_sec - time1.tv_sec;
    }
    NSLog(@"total time %ld", totalCatalogTime);

    gettimeofday(&time1, NULL);
    delete pmap;
    gettimeofday(&time2, NULL);
    NSLog(@"delete pmap time %ld", time2.tv_sec - time1.tv_sec);


    gettimeofday(&time1, NULL);

    [self readPopulations:stacksDocument];

    LocusMO *bLocusMO = [loci.allObjects objectAtIndex:0];
    NSLog(@"pre locus %@ datums %ld", bLocusMO.locusId, bLocusMO.datums.count);

    stacksDocument.loci = loci;
    LocusMO *cLocusMO = [stacksDocument.loci.allObjects objectAtIndex:0];
    NSLog(@"post locus %@ datums %ld", cLocusMO.locusId, cLocusMO.datums.count);


    gettimeofday(&time2, NULL);
    NSLog(@"create stacks document time %ld", time2.tv_sec - time1.tv_sec);


//    gettimeofday(&time1, NULL);
//    NSMutableArray *populations = [stacksDocument findPopulations];
//    NSLog(@"popps %ld", populations.count);
//    gettimeofday(&time2, NULL);
//    NSLog(@"find population time %ld", time2.tv_sec - time1.tv_sec);

    NSLog(@"loading stacks");
    gettimeofday(&time1, NULL);
    [self loadStacks:stacksDocument];
    gettimeofday(&time2, NULL);
    NSLog(@"finished loading stacks time %ld", time2.tv_sec - time1.tv_sec);

    NSError *error;
    BOOL success = [stacksDocument.managedObjectContext save:&error];
    NSLog(@"saved %d", success);

    return stacksDocument;
}

- (void)loadStacks:(StacksDocument *)document {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *path = document.path;
    // TODO: if male . . .male.tags.tsv / female.tags.tsv
    // TODO: or *|!male|!female_N.tags.tsv

    // 1- get all files i the directory named *.tag.tsv"
    NSError *error;
    NSArray *files = [fileManager contentsOfDirectoryAtPath:path error:&error];
    NSLog(@"# of files for directory %ld", files.count);

    // 2 - for each file, read the .tags file
    for (NSString *filePath in files) {
        if ([filePath hasSuffix:@".tags.tsv"] && ![filePath hasPrefix:@"batch"]) {
            [self loadTagFile:document fromFile:filePath];
        }
        else {
            NSLog(@"not loading tag file %@", filePath);
        }
    }


}

- (void)loadTagFile:(StacksDocument *)document fromFile:(NSString *)tagFileName {
    NSLog(@"Loading tag file %@", tagFileName);

    NSUInteger fileNameLength = tagFileName.length;
    NSString *sampleName = [tagFileName substringToIndex:fileNameLength - 9];
    NSLog(@"sampleName %@", sampleName);
    // sampleName . . . from lsat index of "/" . . . to just before ".tags.tsv"

    NSManagedObjectContext *moc = document.managedObjectContext;
    SampleMO* sampleMO = [sampleRepository getSampleForName:sampleName andContext:document.managedObjectContext andError:nil];

    struct timeval time1, time2;
    gettimeofday(&time1, NULL);

    NSError *error2 = nil;
    NSString *absoluteFileName = [document.path stringByAppendingFormat:@"/%@", tagFileName];
    NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error2] componentsSeparatedByString:@"\n"];
    gettimeofday(&time2, NULL);
//    NSLog(@"load file data and split %ld", (time2.tv_sec - time1.tv_sec));

    if (error2) {
        NSLog(@"error loading file [%@]: %@", tagFileName, error2);
    }

//    NSMutableArray *stackEntries = [[NSMutableArray alloc] init];

    NSString *line;
    NSUInteger row = 1;
    gettimeofday(&time1, NULL);
    NSInteger locusId = -1;
    NSInteger newLocusId;
    StackMO *stackMO = nil ;
    DatumMO *datumMO = nil ;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

        if (columns.count > 8) {

            // if the StackMO is found

            newLocusId = [[columns objectAtIndex:2] integerValue];
            if (locusId != newLocusId) {

                // save old
//                NSError* saveError ;
//                BOOL success = [moc save:&saveError];
//                NSLog(@"saved %d", success);
//                if(saveError!=nil ){
//                    NSLog(@"error saving %@",saveError);
//                }

                locusId = newLocusId;
                // search for the new locus
                datumMO = [datumRepository getDatum:moc locusId:locusId andSampleName:sampleMO.name] ;
                if (datumMO!=nil) {
                    if (datumMO.stack == nil) {
                        stackMO = [NSEntityDescription insertNewObjectForEntityForName:@"Stack" inManagedObjectContext:moc];
                        stackMO.datum = datumMO;
                        datumMO.stack = stackMO;
                    }
                    else {
                        stackMO = datumMO.stack;
                    }
                }
                else {
                    datumMO = nil ;
                    stackMO = nil ;
                }
            }

            if (datumMO != nil) {
                StackEntryMO *stackEntryMO = [stackEntryRepository insertStackEntry:moc
                                     entryId:[NSNumber numberWithInteger:row]
                        relationship:[columns objectAtIndex:6]
                        block:[columns objectAtIndex:7]
                        sequenceId:[columns objectAtIndex:8]
                        sequence:[columns objectAtIndex:9]
                        stack:stackMO
                ];

                if ([stackEntryMO.relationship isEqualToString:@"consensus"]) {
                    stackMO.consensus = stackEntryMO;
                }
                else if ([stackEntryMO.relationship isEqualToString:@"model"]) {
                    stackMO.model = stackEntryMO;
                }
                else {
//                    NSLog(@"adding stack entry to stack %ld vs %@",stackMO.datum.locus.locusId,stackMO.datum.sample.name);
                    [stackMO addStackEntriesObject:stackEntryMO];
                    ++row;
                }
            }


        }
    }
    gettimeofday(&time2, NULL);
//    NSLog(@"parse entries lines %ld produce %ld - %ld", fileData.count, stackMO.stackEntries.count, (time2.tv_sec - time1.tv_sec));

    // random snps
//    NSMutableArray *snps = [[NSMutableArray alloc] init];
//
//    for (SnpView *snp in locus.snps) {
//        [snps addObject:[NSNumber numberWithInt:snp.column]];
//    }

//    stacksView.snps = snps;

//    return stacksView;


    // save old
    NSError *saveError;
    BOOL success = [moc save:&saveError];
    NSLog(@"saved %d", success);
    if (saveError != nil ) {
        NSLog(@"error saving %@", saveError);
    }

//    NSEntityDescription *entityDescription3 = [NSEntityDescription entityForName:@"Stack" inManagedObjectContext:moc];
//    NSFetchRequest *request3 = [[NSFetchRequest alloc] init];
//    [request3 setEntity:entityDescription3];
//    NSError *error3;
//    NSArray *stackArray = [moc executeFetchRequest:request3 error:&error3];
//    for (StackMO *stackMO in stackArray) {
//        NSLog(@"# of entries per stack %ld for sample %@ and loci %@", stackMO.stackEntries.count, stackMO.datum.sample.name, stackMO.datum.locus.locusId);
//    }


}

- (void)addSamplesToDocument:(StacksDocument *)document forSampleIds:(vector<int>)sampleIds andSamples:(map<int, string>)samples {

    if (document.populationLookup == nil || document.populationLookup.count == 0) {
        PopulationMO *populationMO = [populationRepository insertPopulation:document.managedObjectContext withId:[NSNumber numberWithInt:1] andName:@"All"];
        document.populations = [NSSet setWithObjects:populationMO, nil];

        // set each sample to populationMO
        for (int i = 0; i < sampleIds.size(); i++) {
            SampleMO *sampleMO = [sampleRepository insertSample:document.managedObjectContext withId:[NSNumber numberWithInt:sampleIds[i]] andName:[NSString stringWithUTF8String:samples[sampleIds[i]].c_str()]];

            [populationMO addSamplesObject:sampleMO];
        }
    }
    else {
        NSNumberFormatter *f = [[NSNumberFormatter alloc] init];
        // else
        for (int i = 0; i < sampleIds.size(); i++) {
            SampleMO *sampleMO = [NSEntityDescription insertNewObjectForEntityForName:@"Sample" inManagedObjectContext:document.managedObjectContext];
            sampleMO.sampleId = [NSNumber numberWithInt:sampleIds[i]];
            sampleMO.name = [NSString stringWithUTF8String:samples[sampleIds[i]].c_str()];

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

    for (NSString *sampleName in populationIdsSet) {
        NSLog(@"adding sample %@", sampleName);
        NSNumber *myNumber = [f numberFromString:sampleName];
        PopulationMO *populationMO = [NSEntityDescription insertNewObjectForEntityForName:@"Population" inManagedObjectContext:document.managedObjectContext];
        populationMO.populationId = myNumber;
        populationMO.name = sampleName;
    }

    NSEntityDescription *entityDescription2 = [NSEntityDescription
            entityForName:@"Population" inManagedObjectContext:document.managedObjectContext];
    NSFetchRequest *request2 = [[NSFetchRequest alloc] init];
    [request2 setEntity:entityDescription2];

    NSError *error;
    NSArray *populationArray = [document.managedObjectContext executeFetchRequest:request2 error:&error];

    document.populations = [NSSet setWithArray:populationArray];

}

- (void)readPopulations:(StacksDocument *)document {
//    NSMutableDictionary *populationLookup = [[NSMutableDictionary alloc] init];
    NSMutableArray *populations = [[NSMutableArray alloc] init];

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
                NSString *sampleName = [columns objectAtIndex:0]; // the sample name . . . male, female, progeny, etc.
                NSManagedObjectContext *moc = document.managedObjectContext;
                NSEntityDescription *entityDescription = [NSEntityDescription entityForName:@"Datum" inManagedObjectContext:moc];
                NSFetchRequest *request = [[NSFetchRequest alloc] init];
                NSPredicate *predicate = [NSPredicate predicateWithFormat:@"name == %@", sampleName];
                [request setPredicate:predicate];
                [request setEntity:entityDescription];
                NSError *error;
                NSArray *datumArray = [moc executeFetchRequest:request error:&error];

                if (datumArray.count == 1) {
//                    DatumMO *datumMO = [datumArray objectAtIndex:0];
                    NSString *populationName = [columns objectAtIndex:1]; // initially an integer
                    NSEntityDescription *entityDescription1 = [NSEntityDescription entityForName:@"Population" inManagedObjectContext:moc];
                    NSFetchRequest *request1 = [[NSFetchRequest alloc] init];
                    NSPredicate *predicate1 = [NSPredicate predicateWithFormat:@"name == %@", populationName];
                    [request1 setPredicate:predicate1];
                    [request1 setEntity:entityDescription1];
                    NSArray *populationArray = [moc executeFetchRequest:request error:&error];

                    PopulationMO *newPopulationMO;
                    if (populationArray == nil || populationArray.count == 0) {
                        newPopulationMO = [NSEntityDescription insertNewObjectForEntityForName:@"Population" inManagedObjectContext:document.managedObjectContext];
                        newPopulationMO.name = populationName;
                        [populations addObject:newPopulationMO];
                    }

                    // [newPopulationMO addDatumsObject:datumMO];
//                    [document.populations insertValue:<#(id)value#> inPropertyWithKey:<#(NSString *)key#>];

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
    NSManagedObjectContext *moc = [stacksDocument getContextForPath:path];
    stacksDocument.managedObjectContext = moc;
    stacksDocument.path = path;
    if (stacksDocumentCreateError) {
        NSLog(@"error creating stacks document %@", stacksDocumentCreateError);
        return nil;
    }
    return stacksDocument;
}

- (StacksDocument *)getStacksDocumentForPath:(NSString *)path {
    NSError *stacksDocumentCreateError;
    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithType:NSSQLiteStoreType error:&stacksDocumentCreateError];
    NSManagedObjectContext *moc = [stacksDocument getContextForPath:path];
//    stacksDocument.managedObjectContext = moc ;
//    stacksDocument.path = path;
    if (stacksDocumentCreateError) {
        NSLog(@"error creating stacks document %@", stacksDocumentCreateError);
        return nil;
    }
    return stacksDocument;
}
@end

