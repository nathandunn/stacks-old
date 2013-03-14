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

#include <dirent.h>
//#include <stdlib.h>


#import "StackEntry.h"
#import "GenotypeView.h"
#import "LocusView.h"
#import "StacksDocument.h"
#import "StacksLoader.h"
#import "StacksView.h"
#import "DataStubber.h"


#include "LociLoader.hpp"

BOOL build_file_list(char const *string1, id param);

@implementation StacksLoader {


}

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

        // TODO: Get the genotype view
        // from populations.cc 150-205 . . . load the catalog matches and then the
        // I will use the locus ID to look over the PopMap<CLocus> and the catalog map map<int,CLocus *> . . .
        // the table that iterates (e.g., is the tabulate_haplotypes object)

        // TODO: get the Stack View
        // different color of view / lgith  / grey is the 3rd column/ locus . . . have to color SNP according other
        // the actual stacksDocuments will need to be imported directly and stored . . . see parse_tsv . .. but will be using raw stacksDocuments

        NSLog(@"model size %d", (int) modelMap.size());

//        stackDocument = [[NSMutableDictionary alloc] initWithCapacity:modelMap.size()];

        map<int, Locus *>::iterator iter = modelMap.begin();
        DataStubber *dataStubber = [[DataStubber alloc] init];

        int randomness = arc4random_uniform(20);
        NSInteger totalGenotypes = 80 + randomness;

        while (iter != modelMap.end()) {
            NSString *sampleId = [NSString stringWithFormat:@"%d", (*iter).first];

            LocusView *locusView = [[LocusView alloc] initWithId:sampleId];

            // TODO: add locus to dictionary / hashmap instead using sampleID as index

            const char *read = (*iter).second->con;
            NSString *letters = [[NSString alloc] initWithCString:read encoding:NSUTF8StringEncoding];
//            NSLog(@"added read %@",letters);
            locusView.consensus = letters;

            // rest of stacksDocuments comes from gentypes . . .  crapola

//            NSMutableArray *generated = [dataStubber generateSnps];
            locusView.snps = [dataStubber generateSnps];
            locusView.male = [dataStubber generateGenotype];
            locusView.female = [dataStubber generateGenotype];
            locusView.progeny = [dataStubber generateProgeny:(NSInteger) totalGenotypes];
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


- (NSMutableDictionary *)loadGenotypes:(NSString *)path withLoci:(NSMutableDictionary *)loci {

    int batch_id = 1;
    //
    // Load the catalog
    //
//    stringstream catalog_file;
    map<int, CSLocus *> catalog;
//    int res;
//    catalog_file << [path UTF8String] << "batch_" << batch_id << ".catalog";
//    if ((res = load_loci(catalog_file.str(), catalog, false)) == 0) {
//        cerr << "Unable to load the catalog '" << catalog_file.str() << "'\n";
//        return 0;
//    }
//
////    NSLog(@"catalog size %d", (int) catalog.size());


//    PopulationLoader* populationLoader = new PopulationLoader();

    // Load matches to the catalog

    vector<vector<CatMatch *> > catalog_matches;
    map<int, string> samples;
    vector<int> sample_ids;

//    srandom(time(NULL));

    vector<pair<int, string> > files;
//    map<int, pair<int, int> > pop_indexes;
    string in_path = [path UTF8String];

//    if (!populationLoader->build_file_list([path UTF8String] ,files, pop_indexes)){
//        exit(1);
//    }
    uint pos;
    string file;
    struct dirent *direntry;

    DIR *dir = opendir(in_path.c_str());

    if (dir == NULL) {
        cerr << "Unable to open directory '" << in_path << "' for reading.\n";
        exit(1);
    }


    while ((direntry = readdir(dir)) != NULL) {
        cout << "reading directory!!!" << endl;
        file = direntry->d_name;

        if (file == "." || file == "..")
            continue;

        if (file.substr(0, 6) == "batch_")
            continue;

        pos = file.rfind(".tags.tsv");
        if (pos < file.length())
            files.push_back(make_pair(1, file.substr(0, pos)));
    }


    cout << "done reading directory!!" << endl;
    cout << "files.size() " << file.size() << endl;


    for (uint i = 0; i < files.size(); i++) {
        vector<CatMatch *> m;
        load_catalog_matches(in_path + files[i].second, m);

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
            exit(0);
        }
    }

    //
    // Create the population map
    //
    cerr << "Populating observed haplotypes for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";
//    PopMap<CLocus> *pmap = new PopMap<CLocus>(sample_ids.size(), catalog.size());
    PopMap<CSLocus> *pmap = new PopMap<CSLocus>(sample_ids.size(), catalog.size());
    pmap->populate(sample_ids, catalog, catalog_matches);


//    map<int, CSLocus *>::iterator it;
//    map<int, ModRes *>::iterator mit;
//    Datum *d;
//    CSLocus *loc;

    // need to load the genotypes in order to get the markers . . .
    //
    // Load the output from the SNP calling model for each individual at each locus. This
    // model output string looks like this:
    //   OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOEOOOOOOEOOOOOOOOOOOOOOOOOOOOOOOOOOOOOUOOOOUOOOOOO
    // and records model calls for each nucleotide: O (hOmozygous), E (hEterozygous), U (Unknown)
    //
//    for (uint i = 0; i < sample_ids.size(); i++) {
//        map<int, ModRes *> modres;
//        load_model_results(in_path + samples[sample_ids[i]], modres);
//
//        if (modres.size() == 0) {
//            cerr << "Warning: unable to find any model results in file '" << samples[sample_ids[i]] << "', excluding this sample from population analysis.\n";
//            continue;
//        }
//
//        for (it = catalog.begin(); it != catalog.end(); it++) {
//            loc = it->second;
//            d = pmap->datum(loc->id, sample_ids[i]);
//
//            if (d != NULL) {
//                if (modres.count(d->id) == 0) {
//                    cerr << "Fatal error: Unable to find model stacksDocuments for catalog locus " << loc->id
//                            << ", sample ID " << sample_ids[i] << ", sample locus " << d->id
//                            << "; likely IDs were mismatched when running pipeline.\n";
//                    exit(0);
//                }
//                d->len = strlen(modres[d->id]->model);
//                d->model = new char [d->len + 1];
//                strcpy(d->model, modres[d->id]->model);
//            }
//        }
//
//        for (mit = modres.begin(); mit != modres.end(); mit++)
//            delete mit->second;
//        modres.clear();
//    }

//    map<int, pair<int, int> > pop_indexes;

//    bool      kernel_smoothed   = false;
//    uint pop_id, start_index, end_index;
//    map<int, pair<int, int> >::iterator pit;
//    stringstream log;
//    log << "batch_" << batch_id << ".populations.log";
//    string log_path = in_path + log.str();

//    ofstream log_fh(log_path.c_str(), ofstream::out);

//    PopSum<CSLocus> *psum = new PopSum<CSLocus>(pmap->loci_cnt(), pop_indexes.size());
//    psum->initialize(pmap);
//
//    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
//        start_index = pit->second.first;
//        end_index   = pit->second.second;
//        pop_id      = pit->first;
//        cerr << "Generating nucleotide-level summary statistics for population " << pop_id << "\n";
//        psum->add_population(catalog, pmap, pop_id, start_index, end_index, log_fh);
//
////        if (kernel_smoothed && loci_ordered) {
////            cerr << "  Generating kernel-smoothed population statistics";
////            if (bootstrap) cerr << " and bootstrap resampling";
////            cerr << "...\n";
////            kernel_smoothed_popstats(catalog, pmap, psum, pop_id, log_fh);
////        }
//    }

//    cerr << "Tallying loci across populations...";
//    psum->tally(catalog);
//    cerr << "done.\n";

    //
    // Idenitfy polymorphic loci, tabulate haplotypes present.
//    LociLoader *lociLoader = new LociLoader();
//    lociLoader->tabulate_haplotypes(catalog, pmap);


//    vector<vector<CatMatch*>>::iterator catalog_match_iterator = catalog_matches.begin();
//    while(catalog_match_iterator!=catalog_matches.end()){
//        vector<CatMatch*> innerVector= (*catalog_match_iterator);
//        vector<CatMatch*>::iterator innerVector_iterator= innerVector.begin();
//        while(innerVector_iterator!=innerVector.end()){
//            CatMatch* match = (*innerVector_iterator);
//            cout << "match: " << match->sample_id;
////            match->
//        }
//
//        cout << "size " << innerVector.size() << endl;
//        catalog_match_iterator++ ;
//    }

//    cout << "sampeid: " << sample_ids.size();
//    map<string, vector<CSLocus*>> *orderedLoci = &(pmap->ordered_loci);
//    cout << "# of ordered loci: " << orderedLoci->size();
//    map<string, vector<CSLocus*>>::iterator iter2 = orderedLoci->begin();
//    while (iter2 != orderedLoci->end()){
//        cout << "key " << iter2->first << endl ;
//        iter2++;
//    }

//    exit(0);

//    NSLog(@"loci size %ld", [loci count]);
//    for (id key in [loci allKeys]) {
//        NSLog(@"%@ - %@", key, [loci objectForKey:key]);
//
//    }

//    map<int, CSLocus *>::iterator iterator = catalog.begin();
//    while (iterator != catalog.end()) {
//        CSLocus *locus = (*iterator).second;
//        GenotypeView *genotypeView = [[GenotypeView alloc] init];
//        // TODO: set Locus in NSDictionary dictionary / hashmap instead using sampleID as index
//
////        locus->
//
////        [returnArray addObject:genotypeView];
//        NSString *key = [NSString stringWithFormat:@"%d", iterator->first];
//        NSLog(@"key %@", key);
//
//        LocusView *locusView = [loci objectForKey:key];
//        NSString *markerString = [NSString stringWithUTF8String:locus->marker.c_str()];
////        cout << "locus model: "<< locus->model << endl ;
//        cout << "locus values marker[" << locus->marker << "] ann[" << locus->annotation << "] con[" << locus->con << "] " << endl;
//        cout << "f [" << locus->f << "] sample[" << locus->sample_id << "]" << endl;
//        NSLog(@"markerString [%@]", markerString);
//        if (markerString != Nil && markerString.length > 0) {
//
//            NSString *newMarker = [NSString stringWithUTF8String:locus->marker.c_str()];
//            NSLog(@"marker [%@]", newMarker);
//            NSLog(@"marker2 [%@]", newMarker);
//            locusView.marker = newMarker;
//            NSLog(@"annotation %@", [NSString stringWithUTF8String:locus->annotation.c_str()]);
//            NSLog(@"con %@", [NSString stringWithUTF8String:locus->con]);
//            NSLog(@"f %f", locus->f);
//            NSLog(@"depth %d", locus->depth);
//            NSLog(@"sample_id %d", locus->sample_id);
//            NSLog(@"trans_gcnt %d", locus->trans_gcnt);
////            PhyLoc phyLoc = locus->loc;
////            NSLog(@"phyLoc %d", phyLoc.bp);
////            NSLog(@"phyLoc %@", [NSString stringWithUTF8String:phyLoc.chr]);
////            NSLog(@"phyLoc %@", phyLoc.strand);
//
//
////            NSLog(@"model %@", [NSString stringWithUTF8String:locus->model]);
//        }
////        locusView.= [NSString stringWithUTF8String:locus->marker.c_str()];
////        [returnArray addObject:genotypeView];
////        locus->marker;
//        iterator++;
//    }

//    exit(0);
//
// TODO: iterate over the catalog . . . -> OR . . . like write_genomic . . . : 736
// as pmap->ordered_loci.begin() . . . etc.
// Datum in pmap->locus(id) contains genotypes table
// genotype is in datum->obshap . . . (size of 1, 2, typically or more)
//    


    return loci;
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
    NSLog(@"number of files %d",files.size());


    // loci loaded . . . now loading genotype

    // TODO: Get the genotype view
    // from populations.cc 150-205 . . . load the catalog matches and then the
    // I will use the locus ID to look over the PopMap<CLocus> and the catalog map map<int,CLocus *> . . .
    // the table that iterates (e.g., is the tabulate_haplotypes object)

    // TODO: get the Stack View
    // different color of view / lgith  / grey is the 3rd column/ locus . . . have to color SNP according other
    // the actual stacksDocuments will need to be imported directly and stored . . . see parse_tsv . .. but will be using raw stacksDocuments

    NSLog(@"model size %d", (int) catalog.size());


    for (uint i = 0; i < files.size(); i++) {
        vector<CatMatch *> m;
        load_catalog_matches([[path stringByAppendingString:@"/"]UTF8String] + files[i].second, m);

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
//    PopMap<CLocus> *pmap = new PopMap<CLocus>(sample_ids.size(), catalog.size());
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

//        locusView.snps = [dataStubber generateSnps];
//        locusView.male = [dataStubber generateGenotype];
//        locusView.female = [dataStubber generateGenotype];
//        locusView.progeny = [dataStubber generateProgeny:(NSInteger) totalGenotypes];
//        locusView.marker = [dataStubber generateMarker];

        [locusViews setObject:locusView forKey:locusView.locusId];
    }


// TODO: iterate over the catalog . . . -> OR . . . like write_genomic . . . : 736
// as pmap->ordered_loci.begin() . . . etc.
// Datum in pmap->locus(id) contains genotypes table
// genotype is in datum->obshap . . . (size of 1, 2, typically or more)





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
