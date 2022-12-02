//
//  AppController.h
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#import <Foundation/Foundation.h>

class Tree;
class CharMatrix;
@class TreeView;

NS_ASSUME_NONNULL_BEGIN

@interface AppController : NSObject {
    // variables here
    Tree* tree;
    IBOutlet TreeView* treeview;
    double turnoverRate;
    double sharingRate;
    double duration;

    CharMatrix* charMatrix;
}

// getters and setters

@property (readwrite) double turnoverRate;
@property (readwrite) double sharingRate;
@property (readwrite) double duration;


- (IBAction)simulate:(id)sender;
- (Tree*)getTree;

@end

NS_ASSUME_NONNULL_END
