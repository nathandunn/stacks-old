//
// Created by ndunn on 3/18/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#import <Foundation/Foundation.h>

//@class LocusView;


@interface LocusCell : NSTableCellView{
}

@property(assign) IBOutlet NSTextField *locusId;
@property(assign) IBOutlet NSTextField *propertyField;
@property(assign) IBOutlet NSTextField *consensusField;

@end