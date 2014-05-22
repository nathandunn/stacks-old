//
//  DatumMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class SampleMO;
@class ColorGenerator;

@interface DatumMO : NSManagedObject{
    int _tagId ;
    BOOL _fetched ;


}
- (void)fetch;


@property (nonatomic, copy) NSString * name;
@property (nonatomic, copy) NSNumber * sampleId;
@property (nonatomic, assign) int tagId;
@property (nonatomic, assign) int  primitiveTagId;
@property (nonatomic, copy) NSString * alleleData;
@property (nonatomic, copy) NSString * depthData;
@property (nonatomic, copy) NSData * haplotypeData;
@property (nonatomic, copy) NSData *snpData;
@property (nonatomic, copy) NSData *stackData;
@property (nonatomic, copy) NSData *metaData;

@property (nonatomic, retain) ColorGenerator *colorGenerator;

- (NSDictionary *)getHaplotypeOrder:(NSMutableArray *)alleleArray;

@end

@interface DatumMO (CoreDataGeneratedAccessors)



- (NSMutableString *)renderHaplotypeHtml:(NSDictionary *)dictionary;

- (NSMutableString *)renderDepthHtml:(NSDictionary *)dictionary;
- (NSMutableString *)renderNameHtml;

- (NSAttributedString *)renderHaplotypes ;
- (NSAttributedString *)renderDepths;
- (NSDictionary*)generateColorForOrder:(NSUInteger) order;


@end
