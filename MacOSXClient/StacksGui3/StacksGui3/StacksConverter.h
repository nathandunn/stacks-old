//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class StacksView;
@class LocusView;
//@class StacksDocument;


@interface StacksConverter : NSObject

- (StacksView *)loadStacksView:(NSString *)filename atPath:(NSString *)path forTag:(NSInteger)tag locus:(LocusView *)locus;
- (NSSet *)loadLociAndGenotypes:(NSString *)path;
- (NSMutableDictionary *)loadPopulation:(NSString *)path;
@end


