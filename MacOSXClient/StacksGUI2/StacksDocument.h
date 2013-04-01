//
//  StacksDocument.h
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "TreeProtocol.h"

@class LocusView;

@interface StacksDocument : NSObject<TreeProtocol>


@property (strong) NSMutableDictionary *locusViews;
@property (strong) NSString *name ;

- (id)initWithLocusView:(NSMutableDictionary*)locusViews;


@end
