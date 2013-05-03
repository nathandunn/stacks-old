//
// Created by NathanDunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

@class StacksView;
@class LocusView;
@class StacksDocument;


@interface StacksConverter : NSObject

@property(nonatomic, strong) NSManagedObjectContext *managedObjectContext;
@property(nonatomic, strong) NSManagedObjectModel *managedObjectModel;

- (id)init ;
- (StacksView *)loadStacksView:(NSString *)filename atPath:(NSString *)path forTag:(NSInteger)tag locus:(LocusView *)locus;
- (StacksDocument *)loadLociAndGenotypes:(NSString *)path;
- (NSMutableDictionary *)loadPopulation:(NSString *)path;

- (StacksDocument *)loadDocument:(StacksDocument *)document;

- (StacksDocument *)createStacksDocumentForPath:(NSString *)path;
@end


