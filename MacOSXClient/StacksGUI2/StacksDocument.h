//
//  StacksDocument.h
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>

@class LocusView;

@interface StacksDocument : NSObject


// TODO: should be an array of locusviews
// this is a wrapper
@property (strong) NSString *locusId ;
@property (strong) LocusView *locusData;

- (id)initWithLocusView:(LocusView *)locusData;


@end
