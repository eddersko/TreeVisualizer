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
    double speciationRate;
    double extinctionRate;
    double duration;
    double sharingRate;
    double alpha;
    double beta;    
    CharMatrix* charMatrix;
}

// getters and setters
@property (readwrite) double speciationRate;
@property (readwrite) double extinctionRate;
@property (readwrite) double duration;
- (IBAction)simulate:(id)sender;
- (Tree*)getTree;

@end

NS_ASSUME_NONNULL_END
