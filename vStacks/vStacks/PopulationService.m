//
// Created by Nathan Dunn on 2/3/14.
// Copyright (c) 2014 University of Oregon. All rights reserved.
//

#import "PopulationService.h"


@implementation PopulationService {

}

+ (PopulationService *)sharedInstance {
    static PopulationService *sharedInstance = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        sharedInstance = [[PopulationService alloc] init];
        // Do any other initialisation stuff here
    });
    return sharedInstance;
}

- (NSString *)validatePopmap:(NSURL *)popmapURL {

    NSFileManager *fileManager = [NSFileManager defaultManager];

    BOOL exists = [fileManager fileExistsAtPath:[popmapURL path]];
    if (!exists) {
        return [NSString stringWithFormat:@"File does not exist at path %@", popmapURL];
    }

    NSError *error;
    @autoreleasepool {
        NSArray *fileData = [[NSString stringWithContentsOfFile:[popmapURL path] encoding:NSUTF8StringEncoding error:&error] componentsSeparatedByString:@"\n"];
        if(error!=nil){
            return [NSString stringWithFormat:@"Error reading file %@: %@", popmapURL, error];
        }
        for (NSString *line in fileData) {
            NSArray *columns = [line componentsSeparatedByString:@"\t"];
            // line.length is because of newline at the end
            if ( line.length > 0 && columns.count != 2) {
                return [NSString stringWithFormat:@"Bad column count %ld for file %@.  Must be two columns separated by a tab.", columns.count, popmapURL];
            }
        }
    }

    return nil;
}
@end