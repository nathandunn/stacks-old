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
// #import "GenotypeView.h"
#import "LocusView.h"
#import "StacksDocument.h"
#import "StacksConverter.h"
#import "StacksView.h"
// #import "DataStubber.h"


// #include "LociLoader.hpp"
#import "GenotypeEntry.h"
#import "SnpView.h"
#import "SnpMO.h"
#import "DatumMO.h"
#import "HaplotypeMO.h"
#import "DepthMO.h"
#import "PopulationMO.h"
#import "SampleMO.h"
//#import "StackEntry.h"


#include <sys/time.h>

@implementation StacksConverter {


}

- (id)init {
    self = [super init];
    if (self) {
        // nothing write now
    }
    return self;
}

- (StacksView *)loadStacksView:(NSString *)filename atPath:(NSString *)path forTag:(NSInteger)tag locus:(LocusView *)locus {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    // TODO: if male . . .male.tags.tsv / female.tags.tsv
    // TODO: or *|!male|!female_N.tags.tsv

    NSString *absoluteFileName = [NSString stringWithFormat:@"%@/%@.tags.tsv", path, filename];

    BOOL existsAtPath = [fileManager fileExistsAtPath:absoluteFileName];
    NSLog(@"absolute file name %@ exists %d", absoluteFileName, existsAtPath);
    if (!existsAtPath) {
        return nil;
    }
    StacksView *stacksView = [[StacksView alloc] init];


    struct timeval time1, time2;
    gettimeofday(&time1, NULL);

    NSError *error = nil;
    NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error] componentsSeparatedByString:@"\n"];
    gettimeofday(&time2, NULL);
    NSLog(@"load file data and split %ld", (time2.tv_sec - time1.tv_sec));

    if (error) {
        NSLog(@"error loading file [%@]: %@", absoluteFileName, error);
    }

    NSMutableArray *stackEntries = [[NSMutableArray alloc] init];

    NSString *line;
    NSUInteger row = 1;
    gettimeofday(&time1, NULL);
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

        if (columns.count > 8 && [[columns objectAtIndex:2] integerValue] == tag) {
            StackEntryMO *stackEntry = [[StackEntryMO alloc] init];
//            stackEntry.entryId = row;
            stackEntry.entryId = [NSNumber numberWithInteger:row];
            stackEntry.relationship = [columns objectAtIndex:6];
            stackEntry.block = [columns objectAtIndex:7];
            stackEntry.sequenceId = [columns objectAtIndex:8];
            stackEntry.sequence = [columns objectAtIndex:9];

            if ([stackEntry.relationship isEqualToString:@"consensus"]) {
                stacksView.consensus = stackEntry;
            }
            else if ([stackEntry.relationship isEqualToString:@"model"]) {
                stacksView.model = stackEntry;
            }
            else {
                ++row;
                [stackEntries addObject:stackEntry];
            }
        }
    }
    gettimeofday(&time2, NULL);
    NSLog(@"parse entries lines %ld produce %ld - %ld", fileData.count, stackEntries.count, (time2.tv_sec - time1.tv_sec));

    stacksView.stackEntries = stackEntries;
    StackEntryMO *referenceStack = [[StackEntryMO alloc] init];
    referenceStack.entryId = 0;
    stacksView.reference = referenceStack;

    // random snps
    NSMutableArray *snps = [[NSMutableArray alloc] init];

    for (SnpView *snp in locus.snps) {
        [snps addObject:[NSNumber numberWithInt:snp.column]];
    }

    stacksView.snps = snps;

    return stacksView;
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
    while (catalogIterator != catalog.end()) {
        LocusMO *locusMO = [NSEntityDescription insertNewObjectForEntityForName:@"Locus" inManagedObjectContext:stacksDocument.managedObjectContext];
        locusMO.locusId = [NSNumber numberWithInt:(*catalogIterator).first];
        const char *read = (*catalogIterator).second->con;
        NSString *letters = [[NSString alloc] initWithCString:read encoding:NSUTF8StringEncoding];
        locusMO.consensus = letters;
        locusMO.marker = [NSString stringWithUTF8String:catalogIterator->second->marker.c_str()];;
        vector<SNP *> snps = catalogIterator->second->snps;
        vector<SNP *>::iterator snpsIterator = snps.begin();

        NSMutableSet *snpsSet = [[NSMutableSet alloc] init];
        for (; snpsIterator != snps.end(); ++snpsIterator) {
            SnpMO *snpMO = [NSEntityDescription insertNewObjectForEntityForName:@"Snp" inManagedObjectContext:stacksDocument.managedObjectContext];
            SNP *snp = (*snpsIterator);
            snpMO.column = [NSNumber numberWithInt:snp->col];
            snpMO.lratio = [NSNumber numberWithFloat:snp->lratio];
            snpMO.rank1 = [NSNumber numberWithChar:snp->rank_1];
            snpMO.rank2 = [NSNumber numberWithChar:snp->rank_2];
            snpMO.rank3 = [NSNumber numberWithChar:snp->rank_3];
            snpMO.rank4 = [NSNumber numberWithChar:snp->rank_4];
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
    for (uint i = 0; i < sample_ids.size(); i++) {
        int sampleId = sample_ids[i];
        string sampleString = samples[sampleId];

        gettimeofday(&time1, NULL);
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
                DatumMO *datumMO = nil ;
                for (DatumMO *datumMO1 in locusMO.datums.allObjects) {
                    if ([datumMO.name isEqualToString:datumMO1.name]) {
                        datumMO = datumMO1;
                    }
                }


                if (datumMO == nil) {
//                    NSLog(@"datum NOT found for key %@ and locus %@",key,locusMO.locusId);
                    vector<char *> obshape = datum->obshap;
                    vector<int> depths = datum->depth;
                    int numLetters = obshape.size();
//                    genotypeMO = [[DatumMO alloc] init];
                    DatumMO *newDatumMO = [NSEntityDescription insertNewObjectForEntityForName:@"Datum" inManagedObjectContext:stacksDocument.managedObjectContext];
                    newDatumMO.name = key;
                    newDatumMO.sampleId = [NSNumber numberWithInt:sample_ids[i]];

                    // get catalogs for matches
                    newDatumMO.tagId = [NSNumber numberWithInt:datum->id];

                    locusMO.length = [NSNumber numberWithInt:loc->depth];

                    if (depths.size() == numLetters) {
                        for (int j = 0; j < numLetters; j++) {
                            HaplotypeMO *haplotypeMO = [NSEntityDescription insertNewObjectForEntityForName:@"Haplotype" inManagedObjectContext:stacksDocument.managedObjectContext];
                            haplotypeMO.haplotype = [NSString stringWithUTF8String:obshape[j]];
                            [newDatumMO addHaplotypesObject:haplotypeMO];

                            DepthMO *depthMO = [NSEntityDescription insertNewObjectForEntityForName:@"Depth" inManagedObjectContext:stacksDocument.managedObjectContext];
                            depthMO.depth = [NSNumber numberWithInt:depths[j]];

                            [newDatumMO addDepthsObject:depthMO];
                        }
                        [locusMO addDatumsObject:newDatumMO];
                    }
                    else {
                        NSLog(@"mismatchon %@", [NSString stringWithUTF8String:sampleString.c_str()]);
                    }
                }
                else {
                    NSLog(@"datum %@ FOUND for key %@ and locus %@", datumMO.name, key, locusMO.locusId);
                }


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


    gettimeofday(&time1, NULL);
    NSMutableArray *populations = [stacksDocument findPopulations];
    NSLog(@"popps %ld", populations.count);
    gettimeofday(&time2, NULL);
    NSLog(@"find population time %ld", time2.tv_sec - time1.tv_sec);

    NSError *error;
    BOOL success = [stacksDocument.managedObjectContext save:&error];
    NSLog(@"saved %d", success);

    return stacksDocument;
}

- (void)addSamplesToDocument:(StacksDocument *)document forSampleIds:(vector<int>)sampleIds andSamples:(map<int, string>)samples {

    if (document.populationLookup == nil || document.populationLookup.count == 0) {
        PopulationMO *populationMO = [NSEntityDescription insertNewObjectForEntityForName:@"Population" inManagedObjectContext:document.managedObjectContext];
        populationMO.populationId = [NSNumber numberWithInt:1];
        populationMO.name = @"All";
        document.populations = [NSSet setWithObjects:populationMO, nil];

        // set each sample to populationMO
        for (int i = 0; i < sampleIds.size(); i++) {
            SampleMO *sampleMO = [NSEntityDescription insertNewObjectForEntityForName:@"Sample" inManagedObjectContext:document.managedObjectContext];
            sampleMO.sampleId = [NSNumber numberWithInt:sampleIds[i]];
            sampleMO.name = [NSString stringWithUTF8String:samples[sampleIds[i]].c_str()];
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

            if(populationId!=nil){
                // lets get the population . . can use lookup, but this is usually pretty small
                NSLog(@"tyring to populate popid %@", populationId);
                for (PopulationMO *populationMO in document.populations) {
                    NSNumber* endNumber = [f numberFromString:populationId];
                    NSLog(@"comparing to %@ vs %@", populationMO.populationId,endNumber);
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
    NSMutableSet* populationIdsSet = [[NSMutableSet alloc] init];
    for(NSString *populationId in populationLookup.allValues){
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
                    DatumMO *datumMO = [datumArray objectAtIndex:0];
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

