//
//  GenomeData.h
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface GenomeData : NSObject

@property (strong) NSString *gene;
//@property (assign) NSString *letters;
@property (strong) NSString *letters;

- (id)initWithGene:(NSString*)gene letters:(NSString*)letters;

@end
