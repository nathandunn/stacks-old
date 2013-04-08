//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class StacksView;
@class StacksDocument;


@interface StacksLoader : NSObject

-(StacksDocument*) loadLoci:(NSString *) path;
- (NSMutableDictionary *)loadGenotypes:(NSString *)path withLoci:(NSMutableDictionary *) loci;

- (StacksView *)loadStacksView:(NSString *)filename atPath:(NSString *)path forTag:(NSInteger)tag locus:(LocusView *)locus;

- (StacksDocument *)loadLociAndGenotypes:(NSString *)path;
@end

//int build_file_list(string in_path, vector<pair<int, string> > &files) ;

