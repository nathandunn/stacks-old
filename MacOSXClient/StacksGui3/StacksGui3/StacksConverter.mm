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
#import "StacksConverter.h"
#import "PopulationRepository.h"
#import "DepthRepository.h"
#import "DatumRepository.h"
#import "HaplotypeRepository.h"
#import "LocusRepository.h"
#import "DatumMO.h"
#import "HaplotypeMO.h"
#import "DepthMO.h"
#import "PopulationMO.h"
#import "SampleMO.h"
#import "SnpRepository.h"
#import "StackEntryRepository.h"
#import "DatumSnpMO.h"
#import "AlleleRepository.h"


#include <sys/time.h>

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
@synthesize datumRepository;
@synthesize depthRepository;
@synthesize haplotypeRepository;
@synthesize locusRepository;
@synthesize populationRepository;
@synthesize sampleRepository;
@synthesize snpRepository;
@synthesize stackEntryRepository;
@synthesize alleleRepository;


// lookups
@synthesize lociDictionary;


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
        alleleRepository= [[AlleleRepository alloc] init];


        lociDictionary = [[NSMutableDictionary alloc] init];
    }
    return self;
}

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

            LocusSnpMO *snpMO = [snpRepository insertLocusSnp:moc
                                                       column:[NSNumber numberWithInt:snp->col]
                                                       lratio:[NSNumber numberWithFloat:snp->lratio]
                                                        rank1:[NSNumber numberWithInt:snp->rank_1]
                                                        rank2:[NSNumber numberWithInt:snp->rank_2]
                                                        rank3:[NSNumber numberWithInt:snp->rank_3]
                                                        rank4:[NSNumber numberWithInt:snp->rank_4]
                                                        locus:locusMO
            ];
            [locusMO addSnpsObject:snpMO];
        }

        [loci addObject:locusMO];
        [lociDictionary setObject:locusMO forKey:locusMO.locusId];
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
            // TODO: use a lookup here to speed up
            for (LocusMO *aLocus in locusArray) {
                NSNumber *lookupKey = [NSNumber numberWithInteger:[[NSString stringWithFormat:@"%ld", it->first] integerValue]];
                if ([lookupKey isEqualToNumber:aLocus.locusId]) {
                    locusMO = aLocus;
                }
            }

            if (datum != NULL && locusMO != nil) {
                NSString *key = [NSString stringWithUTF8String:sampleString.c_str()];

                SampleMO *sampleMO = [sampleRepository getSampleForName:key andContext:moc andError:nil];

                vector<char *> obshape = datum->obshap;
                vector<int> depths = datum->depth;
                int numLetters = obshape.size();
                // TODO: should be entering this for the locus as well?
                DatumMO *newDatumMO = [datumRepository insertDatum:moc name:key sampleId:[NSNumber numberWithInt:sample_ids[i]] sample:sampleMO locus:locusMO];

                // get catalogs for matches
                // TODO: is not the locus the same thing as the id?  can I use the loc->id here?
                newDatumMO.tagId = [NSNumber numberWithInt:datum->id];
                locusMO.length = [NSNumber numberWithInt:loc->depth];

                if (depths.size() == numLetters) {
                    for (int j = 0; j < numLetters; j++) {
                        HaplotypeMO *haplotypeMO = [haplotypeRepository insertHaplotype:moc haplotype:[NSString stringWithUTF8String:obshape[j]]];
                        [newDatumMO addHaplotypesObject:haplotypeMO];

                        DepthMO *depthMO = [depthRepository insertDepth:moc depth:[NSNumber numberWithInt:depths[j]]];
                        [newDatumMO addDepthsObject:depthMO];
                    }
                    [locusMO addDatumsObject:newDatumMO];
                }
                else {
                    NSLog(@"mismatchon %@", [NSString stringWithUTF8String:sampleString.c_str()]);
                }
            }
        }
        gettimeofday(&time2, NULL);
        NSLog(@"iterating sample %d - time %ld", sample_ids[i], (time2.tv_sec - time1.tv_sec));
        totalCatalogTime += time2.tv_sec - time1.tv_sec;


    }

    NSError *innerError = nil ;
    stacksDocument.loci = loci;
    [stacksDocument.managedObjectContext save:&innerError];
    if (innerError != nil) {
        NSLog(@"error doing inner save: ", innerError);
        return nil;
    }

    NSLog(@"total time %ld", totalCatalogTime);

    gettimeofday(&time1, NULL);
    delete pmap;
    gettimeofday(&time2, NULL);
    NSLog(@"delete pmap time %ld", time2.tv_sec - time1.tv_sec);


    gettimeofday(&time1, NULL);
    // TODO: I don't think this does anything hear
    [self readPopulations:stacksDocument];

    LocusMO *bLocusMO = [loci.allObjects objectAtIndex:0];
    NSLog(@"pre locus %@ datums %ld", bLocusMO.locusId, bLocusMO.datums.count);
//
    stacksDocument.loci = loci;
    LocusMO *cLocusMO = [stacksDocument.loci.allObjects objectAtIndex:0];
    NSLog(@"post locus %@ datums %ld", cLocusMO.locusId, cLocusMO.datums.count);

    gettimeofday(&time2, NULL);
    NSLog(@"create stacks document time %ld", time2.tv_sec - time1.tv_sec);


    NSLog(@"loading stack entries");
    gettimeofday(&time1, NULL);
    [self loadStacksEntriesFromTagFile:stacksDocument];
    gettimeofday(&time2, NULL);
    NSLog(@"finished loading stacks entries time %ld", time2.tv_sec - time1.tv_sec);

    NSLog(@"loading snps onto datums");
    gettimeofday(&time1, NULL);
    [self loadSnpsOntoDatum:stacksDocument];
    gettimeofday(&time2, NULL);
    NSLog(@"finished loading snps onto datum time %ld", time2.tv_sec - time1.tv_sec);

    NSLog(@"loading alleles onto datums");
    gettimeofday(&time1, NULL);
    [self loadAllelesOntoDatum:stacksDocument];
    gettimeofday(&time2, NULL);
    NSLog(@"finished loading alleles onto datum time %ld", time2.tv_sec - time1.tv_sec);


    NSError *error;
    BOOL success = [stacksDocument.managedObjectContext save:&error];
    NSLog(@"saved %d", success);

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
            NSLog(@"not loading tag file %@", filePath);
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
            NSLog(@"not loading alleles file %@", filePath);
        }
    }
}

- (void)loadAlleleFileForDatum:(StacksDocument *)document fromFile:(NSString *)alleleFileName {
    NSLog(@"Loading snps file %@", alleleFileName);

    NSUInteger fileNameLength = alleleFileName.length;
    NSString *sampleName = [alleleFileName substringToIndex:fileNameLength - 9];
    NSLog(@"sampleName %@", sampleName);
    // sampleName . . . from lsat index of "/" . . . to just before ".tags.tsv"

    NSManagedObjectContext *moc = document.managedObjectContext;
    SampleMO *sampleMO = [sampleRepository getSampleForName:sampleName andContext:document.managedObjectContext andError:nil];

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
    char allele  ;
    int depth ;
    float ratio;
//    LocusMO *locusMO = nil ;
    NSNumberFormatter *numberFormatter = [[NSNumberFormatter alloc] init];
    numberFormatter.numberStyle = NSNumberFormatterNoStyle;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

        if (columns.count > 4) {

            // if the StackMO is found
//            sampleId = [[columns objectAtIndex:1] integerValue];
            newLocusId = [[columns objectAtIndex:2] integerValue];
            allele = [[columns objectAtIndex:3] charValue];
            ratio = [[columns objectAtIndex:4] floatValue];
            depth = [[columns objectAtIndex:5] intValue];


            if (locusId != newLocusId) {
                locusId = newLocusId;
//                locusMO = [locusRepository getLocus:moc forId:locusId];
                // search for the new locus
                // TODO: get from in-memory lookup?
                datumMO = [datumRepository getDatum:moc locusId:locusId andSampleName:sampleMO.name];
            }

            if (datumMO != nil) {
                DatumAlleleMO* datumAlleleMO = [alleleRepository insertDatumAllele:moc
                                                                ratio:[NSNumber numberWithFloat:ratio]
                                                                 depth:[NSNumber numberWithInt:depth]
                                                                 allele:[NSNumber numberWithInt:allele]
                                                                 datum:datumMO
                ];
//                [datumMO addSnpsObject:datumSnpMO];
                NSLog(@"inserted allele at %@ for sample %@ and locus %@",datumAlleleMO.allele,datumMO.sample.name,datumMO.locus.locusId);
            }
        }
    }
    gettimeofday(&time2, NULL);
    NSLog(@"parse entries lines %ld produce %ld - %ld", fileData.count, datumMO.stackEntries.count, (time2.tv_sec - time1.tv_sec));


    // save old
    NSError *saveError;
    BOOL success = [moc save:&saveError];
    NSLog(@"saved %d", success);
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

    NSManagedObjectContext *moc = document.managedObjectContext;
    SampleMO *sampleMO = [sampleRepository getSampleForName:sampleName andContext:document.managedObjectContext andError:nil];

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
    numberFormatter.numberStyle = NSNumberFormatterNoStyle;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

        if (columns.count > 6) {

            // if the StackMO is found
//            sampleId = [[columns objectAtIndex:1] integerValue];
            newLocusId = [[columns objectAtIndex:2] integerValue];
            column = [[columns objectAtIndex:3] integerValue];
            lratio = [[columns objectAtIndex:4] floatValue];


            if (locusId != newLocusId) {
                locusId = newLocusId;
//                locusMO = [locusRepository getLocus:moc forId:locusId];
                // search for the new locus
                // TODO: get from in-memory lookup?
                datumMO = [datumRepository getDatum:moc locusId:locusId andSampleName:sampleMO.name];
            }

            if (datumMO != nil) {
                DatumSnpMO* datumSnpMO = [snpRepository insertDatumSnp:moc
                                       column:[NSNumber numberWithInteger:column]
                                       lratio:[NSNumber numberWithFloat:lratio]
                                        rank1:[numberFormatter numberFromString:[columns objectAtIndex:5]]
                                        rank2:[numberFormatter numberFromString:[columns objectAtIndex:6]]
                                        rank3:[numberFormatter numberFromString:[columns objectAtIndex:7]]
                                        rank4:[numberFormatter numberFromString:[columns objectAtIndex:8]]
                                        datum:datumMO
                ];
                [datumMO addSnpsObject:datumSnpMO];
//                NSLog(@"inserted snp at %@ for sample %@ and locus %@",datumSnpMO.column,datumMO.sample.name,datumMO.locus.locusId);
            }
        }
    }
    gettimeofday(&time2, NULL);
    NSLog(@"parse entries lines %ld produce %ld - %ld", fileData.count, datumMO.stackEntries.count, (time2.tv_sec - time1.tv_sec));


    // save old
    NSError *saveError;
    BOOL success = [moc save:&saveError];
    NSLog(@"saved %d", success);
    if (saveError != nil ) {
        NSLog(@"error saving %@", saveError);
    }

}

- (void)loadStacksEntriesFromTagFile:(StacksDocument *)document {
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
    SampleMO *sampleMO = [sampleRepository getSampleForName:sampleName andContext:document.managedObjectContext andError:nil];

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

    NSString *line;
    NSUInteger row = 1;
    gettimeofday(&time1, NULL);
    NSInteger locusId = -1;
    NSInteger newLocusId;
    DatumMO *datumMO = nil ;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

        if (columns.count > 8) {

            // if the StackMO is found
            newLocusId = [[columns objectAtIndex:2] integerValue];
            if (locusId != newLocusId) {

                locusId = newLocusId;
                // search for the new locus
                // TODO: get from in-memory lookup?
                datumMO = [datumRepository getDatum:moc locusId:locusId andSampleName:sampleMO.name];
            }

            if (datumMO != nil) {
                NSString *relationship = [columns objectAtIndex:6];

                if ([relationship isEqualToString:@"consensus"]) {
                    datumMO.consensus = [stackEntryRepository insertConsensusStackEntry:moc
                                                                                entryId:[NSNumber numberWithInteger:row]
                                                                                  block:[columns objectAtIndex:7]
                                                                             sequenceId:[columns objectAtIndex:8]
                                                                               sequence:[columns objectAtIndex:9]
                                                                                  datum:datumMO
                    ];
                }
                else if ([relationship isEqualToString:@"model"]) {
                    datumMO.model = [stackEntryRepository insertModelStackEntry:moc
                                                                        entryId:[NSNumber numberWithInteger:row]
                                                                          block:[columns objectAtIndex:7]
                                                                     sequenceId:[columns objectAtIndex:8]
                                                                       sequence:[columns objectAtIndex:9]
                                                                          datum:datumMO
                    ];
                }
                else {
                    StackEntryMO *stackEntryMO = [stackEntryRepository insertStackEntry:moc
                                                                                entryId:[NSNumber numberWithInteger:row]
                                                                           relationship:[columns objectAtIndex:6]
                                                                                  block:[columns objectAtIndex:7]
                                                                             sequenceId:[columns objectAtIndex:8]
                                                                               sequence:[columns objectAtIndex:9]
                                                                                  datum:datumMO
                    ];
//                    [datumMO addStackEntriesObject:stackEntryMO];
                    ++row;
                }
            }
        }
    }
    gettimeofday(&time2, NULL);
    NSLog(@"parse entries lines %ld produce %ld - %ld", fileData.count, datumMO.stackEntries.count, (time2.tv_sec - time1.tv_sec));


    // save old
    NSError *saveError;
    BOOL success = [moc save:&saveError];
    NSLog(@"saved %d", success);
    if (saveError != nil ) {
        NSLog(@"error saving %@", saveError);
    }
}

- (void)addSamplesToDocument:(StacksDocument *)document forSampleIds:(vector<int>)sampleIds andSamples:(map<int, string>)samples {

    if (document.populationLookup == nil || document.populationLookup.count == 0) {
        PopulationMO *populationMO = [populationRepository insertPopulation:document.managedObjectContext id:[NSNumber numberWithInt:1] name:@"All"];
        document.populations = [NSSet setWithObjects:populationMO, nil];

        // set each sample to populationMO
        for (int i = 0; i < sampleIds.size(); i++) {
            SampleMO *sampleMO = [sampleRepository insertSample:document.managedObjectContext id:[NSNumber numberWithInt:sampleIds[i]] name:[NSString stringWithUTF8String:samples[sampleIds[i]].c_str()]];

            [populationMO addSamplesObject:sampleMO];
        }
    }
    else {
        NSNumberFormatter *f = [[NSNumberFormatter alloc] init];
        // else
        for (int i = 0; i < sampleIds.size(); i++) {
            SampleMO *sampleMO = [sampleRepository insertSample:document.managedObjectContext
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
        [populationRepository insertPopulation:document.managedObjectContext id:myNumber name:popId];
    }

    NSArray *populationArray = [populationRepository getAllPopulations:document.managedObjectContext];

    document.populations = [NSSet setWithArray:populationArray];

}

- (void)readPopulations:(StacksDocument *)document {
    NSMutableArray *populations = [[NSMutableArray alloc] init];
    NSManagedObjectContext *moc = document.managedObjectContext;

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
                PopulationMO *newPopulationMO = [populationRepository getPopulation:moc name:populationName];

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
    NSManagedObjectContext *moc = [stacksDocument getContextForPath:path];
    stacksDocument.managedObjectContext = moc;
    stacksDocument.path = path;
    if (stacksDocumentCreateError) {
        NSLog(@"error creating stacks document %@", stacksDocumentCreateError);
        return nil;
    }
    return stacksDocument;
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

