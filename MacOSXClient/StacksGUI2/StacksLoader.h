//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class StacksView;


@interface StacksLoader : NSObject

-(NSMutableDictionary *) loadLoci:(NSString *) path;
- (NSMutableDictionary *)loadGenotypes:(NSString *)path withLoci:(NSMutableDictionary *) loci;
//-(NSMutableArray *) loadGenotypes:(NSString *) path;

- (StacksView *)loadStacksView:(NSString *)string atPath:() path;
@end

