//
//  LocusView.m
//  StacksGUI
//
//  Created by Nathan Dunn on 2/11/13.
//  Copyright (c) 2013 University of Oregon. All rights reserved.
//

#import "LocusView.h"

@implementation LocusView
@synthesize depth ;
@synthesize genotypes;


- (id)initWithId:(NSString *)locusId {
    if ((self = [super init])) {
        self.locusId = locusId;
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

// TODO: fix parental functions
// to identify the # of parents: look at "identify_parents" in genotypes.cc .
- (NSInteger)matchingParents {
    NSInteger count = 0;
//    if (_male) {
//        ++count;
//    }
//    if (_female) {
//        ++count;
//    }
    return count;
}

//- (NSInteger)genotypes {
//    return [_progeny count];
//}

- (BOOL)hasMale {
//    if (_male !=nil ){
//        return (_male.superScript!= nil || _male.subScript !=nil);
//    }
    return FALSE;
}

- (BOOL)hasFemale {
//    if (_female !=nil ){
//        return (_female.superScript!= nil || _female.subScript !=nil);
//    }
    return FALSE;
}

- (NSInteger)count {
    return self.genotypeCount;
//    return self.matchingParents + self.genotypes;
}
@end
