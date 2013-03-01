//
// Created by ndunn on 2/28/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>


@interface GenotypeEntry : NSObject

@property NSInteger  *entryId ;
@property NSString *type;
@property NSString *superScript;
@property NSString *subScript;

// link info
@property NSInteger  *sampleId;
@property NSInteger  *tagId;

@end