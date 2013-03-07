//
//  LocusView.m
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import "LocusView.h"
#import "GenotypeView.h"
#import "StacksView.h"
#import "GenotypeEntry.h"

@implementation LocusView

- (id)initWithId:(NSString *)locusId {
    if ((self = [super init])) {
        self.locusId = locusId;
    }
    return self;
}

- (id)initWithId:(NSString *)locusId consensus:(NSString *)consensus {
    if ((self = [super init])) {
        self.locusId = locusId;
        self.consensus = consensus;
    }
    return self;
}

//@synthesize locusId;
//@synthesize consensus;
//@synthesize marker;
//@synthesize genotypes;
//@synthesize progeny;
//@synthesize ratio;
//@synthesize snps;

- (NSInteger)matchingParents {
    NSInteger count = 0;
    if (_male) {
        ++count;
    }
    if (_female) {
        ++count;
    }
    return count;
}

- (NSInteger)genotypes {
    return [_progeny count];
}

- (BOOL)hasMale {
    if (_male !=nil ){
        return (_male.superScript!= nil || _male.subScript !=nil);
    }
    return FALSE;
}

- (BOOL)hasFemale {
    if (_female !=nil ){
        return (_female.superScript!= nil || _female.subScript !=nil);
    }
    return FALSE;
}

- (NSInteger)count {
 return self.matchingParents + self.genotypes;
}
@end
