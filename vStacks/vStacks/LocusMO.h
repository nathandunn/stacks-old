//
//  LocusMO.h
//  StacksGui3
//
//  Created by Nathan Dunn on 5/15/13.
//  Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

@class ColorGenerator;

@interface LocusMO : NSManagedObject

@property (nonatomic, copy) NSNumber * basePairs;
@property (nonatomic, copy) NSString * chromosome;
@property (nonatomic, copy) NSString * strand;
@property (nonatomic, copy) NSString * type;
@property (nonatomic, copy) NSString * consensus;
@property (nonatomic, copy) NSNumber * length;
@property (nonatomic, copy) NSNumber * locusId;
@property (nonatomic, copy) NSNumber * parentCount;
@property (nonatomic, copy) NSNumber * progenyCount;
@property (nonatomic, copy) NSString * marker;
@property (nonatomic, copy) NSString * ratio;
@property (nonatomic, copy) NSData *metaData;

@property (nonatomic, copy) NSData *alleleData;
@property (nonatomic, copy) NSData *snpData;

// not a managed part, populated by datum rendering . . .
@property (atomic, retain) ColorGenerator *colorGenerator;



@end

@interface LocusMO (CoreDataGeneratedAccessors)

@end
