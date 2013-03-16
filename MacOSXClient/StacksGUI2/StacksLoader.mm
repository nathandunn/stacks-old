//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
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

BOOL build_file_list(char const *string1, id param);

@implementation StacksLoader {


}

/**
* TODO: remove once other is implmemented
*/
- (StacksDocument *)loadLoci:(NSString *)examplePath {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];

    NSMutableDictionary *locusViews = [[NSMutableDictionary alloc] init];
    if (!existsAtPath) {
        NSLog(@"files do not exist %@", examplePath);
        exit(0);
    }
    else {
//        map<int,ModRes*> modelMap ;
        map<int, Locus *> modelMap;
        NSString *exampleFile = [examplePath stringByAppendingString:@"batch_1.catalog"];

        load_loci([exampleFile UTF8String], modelMap, false);

        NSLog(@"model size %d", (int) modelMap.size());

//        stackDocument = [[NSMutableDictionary alloc] initWithCapacity:modelMap.size()];

        map<int, Locus *>::iterator iter = modelMap.begin();
        DataStubber *dataStubber = [[DataStubber alloc] init];

        int randomness = arc4random_uniform(20);
        NSInteger totalGenotypes = 80 + randomness;

        while (iter != modelMap.end()) {
            NSString *sampleId = [NSString stringWithFormat:@"%d", (*iter).first];

            LocusView *locusView = [[LocusView alloc] initWithId:sampleId];

            const char *read = (*iter).second->con;
            NSString *letters = [[NSString alloc] initWithCString:read encoding:NSUTF8StringEncoding];
            locusView.consensus = letters;

            // rest of stacksDocuments comes from gentypes . . .  crapola
            locusView.snps = [dataStubber generateSnps];
            locusView.genotypes = [dataStubber generateGenotypes:(NSInteger) totalGenotypes];
            locusView.marker = [dataStubber generateMarker];


            [locusViews setObject:locusView forKey:locusView.locusId];

//            StacksDocument *doc = [[StacksDocument alloc] initWithLocusView:locusView];
//            [stackDocuments insertValue:doc atIndex:doc inPropertyWithKey:<#(NSString *)key#>]
//            [stackDocuments insertValue:doc inPropertyWithKey:doc.locusId];
//            [stackDocument setObject:doc forKey:doc.locusId];

            ++iter;
        }
    }

    StacksDocument *stackDocument = [[StacksDocument alloc] initWithLocusView:locusViews];

    return stackDocument;

}


- (StacksView *)loadStacksView:(NSString *)stackKey atPath:(NSString *)path {
    NSFileManager *fileManager = [NSFileManager defaultManager];
    // TODO: if male . . .male.tags.tsv / female.tags.tsv
    // TODO: or *|!male|!female_N.tags.tsv
    NSString *fileName = [[NSString alloc] init];
    if ([stackKey isEqualToString:@"male"]) {
        fileName = @"/male.tags.tsv";
    }
    else if ([stackKey isEqualToString:@"female"]) {
        fileName = @"/female.tags.tsv";
    }
    else {
        fileName = [NSString stringWithFormat:@"/progeny_%@.tags.tsv", stackKey];
    }
    NSString *absoluteFileName = [path stringByAppendingString:fileName];

    BOOL existsAtPath = [fileManager fileExistsAtPath:absoluteFileName];
    NSLog(@"absolute file name %@ exists %d", absoluteFileName, existsAtPath);
    if (!existsAtPath) {
        return nil;
    }
    StacksView *stacksView = [[StacksView alloc] init];

    NSError *error = nil ;
    NSArray *fileData = [[NSString stringWithContentsOfFile:absoluteFileName encoding:NSUTF8StringEncoding error:&error] componentsSeparatedByString:@"\n"];

    if (error) {
        NSLog(@"error loading file [%@]: %@", absoluteFileName, error);
    }

    NSMutableArray *stackEntries = [[NSMutableArray alloc] init];

    NSString *line;
    NSUInteger row = 1;
    for (line in fileData) {
        NSArray *columns = [line componentsSeparatedByString:@"\t"];

        if ([[columns objectAtIndex:0] isEqualToString:@"0"]) {
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

//        for(column in columns){
////            NSLog(@"entry %@",column);
//        }
        }


    }

    stacksView.stackEntries = stackEntries;
    StackEntry *referenceStack = [[StackEntry alloc] init];
    referenceStack.sequence = @"1203130120321302103210321203";
    referenceStack.entryId = 0;
    stacksView.reference = referenceStack;

    // random snps
    NSMutableArray *snps = [[NSMutableArray alloc] init];
    [snps addObject:[NSNumber numberWithInt:17]];
    [snps addObject:[NSNumber numberWithInt:35]];
    stacksView.snps = snps;


    return stacksView;
}

- (StacksDocument *)loadLociAndGenotypes:(NSString *)path {
    [self checkFile:path];
//    map<int, Locus *> modelMap;
    map<int, CSLocus *> catalog;
    NSString *exampleFile = [path stringByAppendingString:@"batch_1.catalog"];

//        load_model_results([exampleFile UTF8String], modelMap);
    load_loci([exampleFile UTF8String], catalog, false);

    /**
    * START: loading genotypes + extra info
*/

    vector<vector<CatMatch *> > catalog_matches;
    map<int, string> samples;
    vector<int> sample_ids;

    vector<pair<int, string>> files = [self buildFileList:path];
//    vector<pair<int, string>>::iterator fileNameIterator = files.begin();
//    while (fileNameIterator != files.end()) {
//        string fileName = (*fileNameIterator).second;
//        NSLog(@"filestring %@", [NSString stringWithUTF8String:fileName.c_str()]);
////        NSLog(@"filestring %@", fileName);
//        ++fileNameIterator;
//    }
    NSLog(@"number of files %d", files.size());


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
        NSString *sampleId = [NSString stringWithFormat:@"%d", (*catalogIterator).first];

        LocusView *locusView = [[LocusView alloc] initWithId:sampleId];
        const char *read = (*catalogIterator).second->con;
        NSString *letters = [[NSString alloc] initWithCString:read encoding:NSUTF8StringEncoding];
//            NSLog(@"added read %@",letters);
        locusView.consensus = letters;
        locusView.marker = [NSString stringWithUTF8String:catalogIterator->second->marker.c_str()];;
        locusView.genotypeCount = catalogIterator->second->gcnt;
        vector<SNP *> snps = catalogIterator->second->snps;
        vector<SNP *>::iterator snpsIterator = snps.begin();

        NSMutableArray *snpsArray = [[NSMutableArray alloc] initWithCapacity:snps.size()];
        for (; snpsIterator != snps.end(); ++snpsIterator) {
            [snpsArray addObject:[NSValue valueWithPointer:(*snpsIterator)]];
        }

//        locusView.male = [dataStubber generateGenotype];
//        locusView.female = [dataStubber generateGenotype];
//        locusView.progeny = [dataStubber generateProgeny:(NSInteger) totalGenotypes];

        [locusViews setObject:locusView forKey:locusView.locusId];
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

            LocusView *locusView = [locusViews objectForKey:[NSString stringWithFormat:@"%d", it->first]];
            if (d != NULL && locusView != nil) {
//                NSString *key = [NSString stringWithFormat:@"%d",sample_ids[i]];
                NSString *key = [NSString stringWithUTF8String:sampleString.c_str()];
                NSMutableDictionary *genotypes = locusView.genotypes;
                if(genotypes==nil){
                    genotypes=[[NSMutableDictionary alloc] init];
                }
                GenotypeEntry *genotypeEntry = [genotypes objectForKey:key];
                if (genotypeEntry == nil) {
                    genotypeEntry = [[GenotypeEntry alloc] init];
                }
                genotypeEntry.name = key;
                genotypeEntry.sampleId = [[NSNumber numberWithInt:sample_ids[i]] unsignedIntegerValue];

                locusView.depth = loc->depth;

//                NSLog(@"locus: %d sample %d",it->first,sample_ids[i]) ;
//                NSLog(@"id: %d",d->id);
//                NSLog(@"length: %d",d->len);
//                NSLog(@"tot_depth: %d",d->tot_depth);

//                NSLog(@"objshape size: %d",obshape.size());
                vector<char *> obshape = d->obshap;
                vector<int> depths = d->depth;
                int numLetters = obshape.size();
                if (depths.size() == numLetters) {
                    genotypeEntry.haplotypes = [[NSMutableArray alloc] initWithCapacity:numLetters];
                    for (int j = 0; j < numLetters; j++) {
                        [genotypeEntry.haplotypes addObject:[NSString stringWithUTF8String:obshape[j]]];
                        [genotypeEntry.depths addObject:[NSNumber numberWithInt:depths[j]]];
                    }
                    [genotypes setObject:genotypeEntry forKey:key];
                    locusView.genotypes = genotypes;
                }
                else {
                    NSLog(@"mismatchon %@", [NSString stringWithUTF8String:sampleString.c_str()]);
//                    sampleString);
                }


            }
        }
    }


    StacksDocument *stackDocument = [[StacksDocument alloc] initWithLocusView:locusViews];

    return stackDocument;
}

/**
* returns all of the files that end with .tags.tsv in the directory
*/
- (vector<pair<int, string> >)buildFileList:(NSString *)path {
    vector<pair<int, string> > files;

    NSFileManager *fileManager = [NSFileManager defaultManager];

    NSDirectoryEnumerator *dirEnum = [fileManager enumeratorAtPath:path];
    NSString *file;
    while (file = [dirEnum nextObject]) {

        if ([file hasSuffix:@".tags.tsv"]) {
//            NSLog(@"HAS SUFFIX name %@", file);
            NSUInteger length = [file length] - 9;
            NSString *fileName = [file substringToIndex:length];
            files.push_back(make_pair(1, [fileName UTF8String]));
        }
    }

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

//int build_file_list(string in_path, vector<pair<int, string> > &files) {
//    uint pos;
//    string file;
//    struct dirent *direntry;
//
//    DIR *dir = opendir(in_path.c_str());
//
//    if (dir == NULL) {
//        cerr << "Unable to open directory '" << in_path << "' for reading.\n";
//        exit(1);
//    }
//
//
//    while ((direntry = readdir(dir)) != NULL) {
//        cout << "reading directory!!!" << endl;
//        file = direntry->d_name;
//
//        if (file == "." || file == "..")
//            continue;
//
//        if (file.substr(0, 6) == "batch_")
//            continue;
//
//        pos = file.rfind(".tags.tsv");
//        if (pos < file.length())
//            files.push_back(make_pair(1, file.substr(0, pos)));
//    }
//
//
//    cout << "done reading directory!!" << endl;
//    cout << "files.size() " << file.size() << endl;
//    return files.size();
//}
