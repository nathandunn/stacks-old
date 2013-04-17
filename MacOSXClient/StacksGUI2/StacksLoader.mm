//
// Created by NathanDunn on 2/28/13.
//
//



#include "locus.h"
//#include "stacks.h"
#include "sql_utilities.h"
#include "PopMap.h"
#import "StacksDocument.h"


#include <fstream>

using std::ifstream;
using std::ofstream;

//#include "PopSum.h"

//#include <dirent.h>
//#include <stdlib.h>


#import "StackEntry.h"
#import "GenotypeView.h"
#import "LocusView.h"
#import "StacksDocument.h"
#import "StacksLoader.h"
#import "StacksView.h"
#import "DataStubber.h"


#include "LociLoader.hpp"
#import "GenotypeEntry.h"
#import "SnpView.h"

@implementation StacksLoader {


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

    NSError *error = nil;
    NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error] componentsSeparatedByString:@"\n"];

    if (error) {
        NSLog(@"error loading file [%@]: %@", absoluteFileName, error);
    }

    NSMutableArray *stackEntries = [[NSMutableArray alloc] init];

    NSString *line;
    NSUInteger row = 1;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

//        NSLog(@"columns at index %@",[columns objectAtIndex:2]);
        if (columns.count > 8 && [[columns objectAtIndex:2] integerValue] == tag) {
            StackEntry *stackEntry = [[StackEntry alloc] init];
            stackEntry.entryId = row;
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

    stacksView.stackEntries = stackEntries;
    StackEntry *referenceStack = [[StackEntry alloc] init];
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
    [self checkFile:path];
    map<int, CSLocus *> catalog;
    NSString *catalogFile = [path stringByAppendingString:@"batch_1.catalog"];
    load_loci([catalogFile UTF8String], catalog, false);
    NSMutableDictionary *populationLookup = [self loadPopulation:path];

    /**
    * START: loading genotypes + extra info
*/
    vector<vector<CatMatch *> > catalog_matches;
    map<int, string> samples;
    vector<int> sample_ids;

    vector<pair<int, string>> files = [self buildFileList:path];
    NSLog(@"number of files %ld", files.size());


    // loci loaded . . . now loading genotype
    NSLog(@"model size %d", (int) catalog.size());

    for (uint i = 0; i < files.size(); i++) {
        vector<CatMatch *> m;
        load_catalog_matches([[path stringByAppendingString:@"/"] UTF8String] + files[i].second, m);

        if (m.size() == 0) {
            cerr << "Warning: unable to find any matches in file '" << files[i].second << "', excluding this sample from population analysis.\n";
//            exit(1);
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


    cerr << "Populating observed haplotypes for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";
    PopMap<CSLocus> *pmap = new PopMap<CSLocus>(sample_ids.size(), catalog.size());
    pmap->populate(sample_ids, catalog, catalog_matches);


    NSMutableDictionary *locusViews = [[NSMutableDictionary alloc] init];

    map<int, CSLocus *>::iterator catalogIterator = catalog.begin();
    while (catalogIterator != catalog.end()) {
        LocusView *locusView = [[LocusView alloc] initWithId:(*catalogIterator).first];
        const char *read = (*catalogIterator).second->con;
        NSString *letters = [[NSString alloc] initWithCString:read encoding:NSUTF8StringEncoding];
        locusView.consensus = letters;
        locusView.marker = [NSString stringWithUTF8String:catalogIterator->second->marker.c_str()];;
        vector<SNP *> snps = catalogIterator->second->snps;
        vector<SNP *>::iterator snpsIterator = snps.begin();

        NSMutableArray *snpsArray = [[NSMutableArray alloc] initWithCapacity:snps.size()];
        for (; snpsIterator != snps.end(); ++snpsIterator) {
            SnpView *snpView = [[SnpView alloc] init];
            SNP *snp = (*snpsIterator);
            snpView.column = snp->col;
            snpView.lratio = snp->lratio;
            snpView.rank1 = snp->rank_1;
            snpView.rank2 = snp->rank_2;
            snpView.rank2 = snp->rank_2;
            snpView.rank4 = snp->rank_4;
            [snpsArray addObject:snpView];
        }

        locusView.snps = snpsArray;

        [locusViews setObject:locusView forKey:[NSString stringWithFormat:@"%ld",locusView.locusId] ];
        ++catalogIterator;
    }


// TODO: iterate over the catalog . . . -> OR . . . like write_genomic . . . : 736
// as pmap->ordered_loci.begin() . . . etc.
// Datum in pmap->locus(id) contains genotypes table
// genotype is in datum->obshap . . . (size of 1, 2, typically or more)

    map<int, CSLocus *>::iterator it;
    Datum *d;
    CSLocus *loc;

    // for each sample process the catalog

    for (uint i = 0; i < sample_ids.size(); i++) {
        int sampleId = sample_ids[i];
        string sampleString = samples[sampleId];

//        NSLog(@"evaluating sample %@", [NSString stringWithUTF8String:sampleString.c_str()]);
        for (it = catalog.begin(); it != catalog.end(); it++) {
            loc = it->second;
            d = pmap->datum(loc->id, sample_ids[i]);

            LocusView *locusView = [locusViews objectForKey:[NSString stringWithFormat:@"%ld", it->first]];
            if (d != NULL && locusView != nil) {
//                NSString *key = [NSString stringWithFormat:@"%d",sample_ids[i]];
                NSString *key = [NSString stringWithUTF8String:sampleString.c_str()];
                NSMutableDictionary *genotypes = locusView.genotypes;
                if (genotypes == nil) {
                    genotypes = [[NSMutableDictionary alloc] init];
                }
                GenotypeEntry *genotypeEntry = [genotypes objectForKey:key];
                if (genotypeEntry == nil) {
                    genotypeEntry = [[GenotypeEntry alloc] init];
                }
                genotypeEntry.name = key;
                genotypeEntry.sampleId = [[NSNumber numberWithInt:sample_ids[i]] unsignedIntegerValue];

                // get catalogs for matches
                genotypeEntry.tagId = d->id;

                locusView.depth = loc->depth;

                vector<char *> obshape = d->obshap;
                vector<int> depths = d->depth;
                int numLetters = obshape.size();
                if (depths.size() == numLetters) {
                    genotypeEntry.haplotypes = [[NSMutableArray alloc] initWithCapacity:numLetters];
                    genotypeEntry.depths = [[NSMutableArray alloc] initWithCapacity:numLetters];
                    for (int j = 0; j < numLetters; j++) {
                        [genotypeEntry.haplotypes addObject:[NSString stringWithUTF8String:obshape[j]]];
                        [genotypeEntry.depths addObject:[NSNumber numberWithInt:depths[j]]];
                    }
//                    for(NSNumber *depth in genotypeEntry.depths){
//                        NSLog(@"depth %@",depth);
//                    }
                    [genotypes setObject:genotypeEntry forKey:key];
                    locusView.genotypes = genotypes;
                }
                else {
                    NSLog(@"mismatchon %@", [NSString stringWithUTF8String:sampleString.c_str()]);
                }
            }
        }
    }

    delete pmap ;


    StacksDocument *stacksDocument = [[StacksDocument alloc] initWithLocusView:locusViews];
    stacksDocument.path = path;
    stacksDocument.populationLookup = populationLookup ;

    NSMutableArray *populations = [stacksDocument findPopulations];
    NSLog(@"popps %ld",populations.count);

    return stacksDocument;
}

- (NSMutableDictionary *)loadPopulation:(NSString *)path {
    NSMutableDictionary *populationLookup = [[NSMutableDictionary alloc] init];

    NSString *popmapFile = [path stringByAppendingString:@"popmap"];
    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL exists = [fileManager fileExistsAtPath:popmapFile];

    if(exists){
        NSArray *fileData = [[NSString stringWithContentsOfFile:popmapFile encoding:NSUTF8StringEncoding error:nil] componentsSeparatedByString:@"\n"];
        NSString *line ;
        for(line in fileData){
            NSArray *columns = [line componentsSeparatedByString:@"\t"];
            if(columns.count==2){
                [populationLookup setObject:[columns objectAtIndex:1] forKey:[columns objectAtIndex:0]];
            }
            else{
                NSLog(@"something wrong with the column count %@",line);
            }
        }
    }
    else{
        NSLog(@"does not exist at %@",popmapFile);
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

@end

