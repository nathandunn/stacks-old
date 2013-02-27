//
//  StacksDocument.h
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>

@class GenomeData;

@interface StacksDocument : NSObject

@property (strong) GenomeData *data;

- (id)initWithGene:(NSString*)gene letters:(NSString*)letters ;


@end
