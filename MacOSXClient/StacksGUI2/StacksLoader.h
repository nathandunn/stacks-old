//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>


@interface StacksLoader : NSObject

-(NSMutableArray *) loadLoci:(NSString *) path;
-(NSMutableArray *) loadGenotypes:(NSString *) path;

@end

