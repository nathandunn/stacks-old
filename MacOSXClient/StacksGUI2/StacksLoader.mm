//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StacksDocument.h"
#import "locus.h"
#import "stacks.h"
#import "sql_utilities.h"
#import "StacksLoader.h"
#import "LocusView.h"
#import "PopMap.h"
//#import "CLocus.hpp"
//#import "PopulationLoader.hpp"
#import "GenotypeView.h"




@implementation StacksLoader {
    

}

- (NSMutableArray *)loadLoci:(NSString *)examplePath{
    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];

    // TODO: should be NSMutableDictionary
    NSMutableArray *geneDocs = [[NSMutableArray alloc] init];
    if(!existsAtPath){
        NSLog(@"files do not exist %@",examplePath);
        StacksDocument *doc1 = [[StacksDocument alloc] initWithMarker:@"sox9" consensusSequence:@"ATATAGATA"];
        StacksDocument *doc2 = [[StacksDocument alloc] initWithMarker:@"sox12" consensusSequence:@"ATATAGAGG"];
        geneDocs = [NSMutableArray arrayWithObjects:doc1,doc2,nil];
    }
    else{
//        map<int,ModRes*> modelMap ;
        map<int,Locus*> modelMap ;
        NSString *exampleFile = [examplePath  stringByAppendingString:@"batch_1.catalog"] ;

//        load_model_results([exampleFile UTF8String], modelMap);
        load_loci([exampleFile UTF8String], modelMap,false);

        // TODO: Get the genotype view
        // from populations.cc 150-205 . . . load the catalog matches and then the
        // I will use the locus ID to look over the PopMap<CLocus> and the catalog map map<int,CLocus *> . . .
        // the table that iterates (e.g., is the tabulate_haplotypes object)

        // TODO: get the Stack View
        // different color of view / lgith  / grey is the 3rd column/ locus . . . have to color SNP according other
        // the actual data will need to be imported directly and stored . . . see parse_tsv . .. but will be using raw data

        NSLog(@"model size %d",(int)modelMap.size());

//        NSArray *dirFiles = [fileManager contentsOfDirectoryAtPath:examplePath error:nil];
//        NSArray *fastaFiles= [dirFiles filteredArrayUsingPredicate:[NSPredicate predicateWithFormat:@"self ENDSWITH '.fa'"]];
//        NSEnumerator *e = [fastaFiles objectEnumerator];
//        id object;
        geneDocs = [NSMutableArray arrayWithCapacity:modelMap.size()];
        map<int,Locus*>::iterator iter = modelMap.begin();

        while(iter!=modelMap.end()){
            NSString *sampleId = [NSString stringWithFormat:@"%d",(*iter).first];
            LocusView *locusView = [[LocusView alloc] initWithId:sampleId ];

            // TODO: add locus to dictionary / hashmap instead using sampleID as index

            const char *read = (*iter).second->con;
            NSString* letters = [[NSString alloc] initWithCString:read encoding: NSUTF8StringEncoding];
//            NSLog(@"added read %@",letters);
            locusView.consensus = letters;
            
            // rest of data comes from gentypes . . .  crapola



            StacksDocument *doc = [[StacksDocument alloc] initWithLocusData:locusView];
            [geneDocs addObject:doc];

            ++iter;
        }
    }
    
    return geneDocs;

}


//- (NSMutableArray *)loadLoci:(NSString *)examplePath{
- (NSMutableArray *)loadGenotypes:(NSString *)path withLoci:(NSMutableArray *) loci{
    
    NSMutableArray *returnArray = [[NSMutableArray alloc] init];

    int batch_id = 1 ;
    //
    // Load the catalog
    //
    stringstream catalog_file;
    map<int, CSLocus *> catalog;
    int res;
    catalog_file << [path UTF8String] << "batch_" << batch_id << ".catalog";
    if ((res = load_loci(catalog_file.str(), catalog, false)) == 0) {
        cerr << "Unable to load the catalog '" << catalog_file.str() << "'\n";
        return 0;
    }

    NSLog(@"catalog size %d",(int) catalog.size());


//    PopulationLoader* populationLoader = new PopulationLoader();

    // Load matches to the catalog
    //
    vector<vector<CatMatch *> > catalog_matches;
    map<int, string>            samples;
    vector<int>                 sample_ids;

    srandom(time(NULL));

    vector<pair<int, string> > files;
//    map<int, pair<int, int> > pop_indexes;
    string in_path ;

//    if (!populationLoader->build_file_list([path UTF8String] ,files, pop_indexes)){
//        exit(1);
//    }


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
            cerr << "Fatal error: sample ID " << m[0]->sample_id << " occurs twice in this data set, likely the pipeline was run incorrectly.\n";
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

    map<int,CSLocus*>::iterator iterator= catalog.begin();

    while(iterator!=catalog.end()){
        CSLocus* locus = (*iterator).second;
        GenotypeView * genotypeView = [[GenotypeView alloc] init];
        // TODO: set Locus in NSDictionary dictionary / hashmap instead using sampleID as index

//        locus->

        [returnArray addObject:genotypeView];
//        locus->marker;
        iterator++;
    }


    return returnArray;
}


@end