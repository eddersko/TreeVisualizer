//
//  TreeView.h
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#import <Cocoa/Cocoa.h>
@class AppController;

NS_ASSUME_NONNULL_BEGIN

@interface TreeView : NSView {
    
    IBOutlet AppController* appController;
    
}

@end

NS_ASSUME_NONNULL_END
