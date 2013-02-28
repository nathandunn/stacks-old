//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import "StacksDocument.h"
#import "stacks.h"
#import "sql_utilities.h"
#import "StacksLoader.h"


@implementation StacksLoader {

}

- (NSMutableArray *)loadLoci:(NSString *)examplePath{
    NSFileManager *fileManager = [NSFileManager defaultManager];
    BOOL existsAtPath = [fileManager fileExistsAtPath:examplePath];
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
            const char *read = (*iter).second->con;
//            NSString* letters = [NSString stringWithCString:read encoding: NSUTF8StringEncoding];
//            [NSString initWithUTF8String:read];

//            NSString* letters = [NSString stringWithUTF8String:read];
//            const char *test = "teststring" ;
            NSString* letters = [[NSString alloc] initWithCString:read encoding: NSUTF8StringEncoding];
            //            NSString* letters = [NSString stringWithFormat:@"%@",lettersString];
//            NSString* letters = [NSString stringWithCharacters:read length: sizeof(read)];
//            NSString* letters = [NSString stringWithCharacters:read length: sizeof(read)];

//            NSString* letters = @"";
//            letters = @"ABCDEFG";
            NSLog(@"added read %@",letters);
            int intValue = (*iter).first;
            NSString *geneValue = [NSString stringWithFormat:@"%d",intValue];
            StacksDocument *doc = [[StacksDocument alloc] initWithMarker:geneValue consensusSequence:letters];
            [geneDocs addObject:doc];

            ++iter;
        }

//        while (object = [e nextObject]) {
//    //        NSLog(object);
//            // do something with object
//        }
    }
    
    return geneDocs;

}


@end