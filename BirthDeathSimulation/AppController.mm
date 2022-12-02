//
//  AppController.m
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#import "AppController.h"
#include "Tree.hpp"
#import "TreeView.h"
#include "CharMatrix.hpp"
#include "RandomVariable.hpp"

@implementation AppController;
@synthesize speciationRate;
@synthesize extinctionRate;
@synthesize duration;

- (IBAction)simulate:(id)sender {
    
    NSLog(@"In simulate function.");
    if (tree != NULL)
        delete tree;
    
    // simulate new tree
    tree = new Tree(speciationRate, extinctionRate, duration);
    
    // tree with ghost lineages
    //tree = new Tree("(((A:0.1,B:0.1):0.1,E:0.1):0.1,(C:0.2,D:0.2):0.1);");
    
    // tree without ghost lineages
    //tree = new Tree("((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);");

    // begin initializing the rate matrix
    int numStates = 2;
    
    //q = array of pointers to doubles with size ns (2)
    double** q = new double*[numStates];
    
    // q[0] = array of doubles with size ns*ns (4)
    q[0] = new double[numStates*numStates];
    
    // assign memory address to q[1]
    for (int i = 1; i < numStates; i++)
        q[i] = q[i-1]+numStates;
    
    // f = array of doubles, equilibrium frequencies?
    double* f = new double[numStates];
    
    // f[0] = 0.5, f[1] = 0.5
    for (int i = 0; i < numStates; i++)
        f[i] = 1.0/numStates;
    
    // q[0][1], q[1][0] = 1.0
    // [] [1.0]
    // [1.0] []
    for (int i = 0; i < numStates; i++) {
        for (int j = i+1; j < numStates; j++) {
            q[i][j] = 1.0;
            q[j][i] = 1.0;
        }
    }
    
    // diagonals is negative sum of row
    // q[0][0], q[1][1] = -1.0
    // [-1.0] [1.0]
    // [1.0] [-1.0]
    for (int i = 0; i < numStates; i++) {
        double sum = 0.0;
        for (int j = 0; j < numStates; j++) {
            if (i != j)
                sum += q[i][j];
        }
        q[i][i] = -sum;
    }
    
    double averageRate = 0.0;
    
    // average rate is (1.0*0.5) + (1.0 * 0.5) = 1.0
    for (int i = 0; i < numStates; i++) {
        averageRate += -q[i][i] * f[i];
    }
    
    // scaleFactor = 1.0
    double scaleFactor = 1.0 / averageRate;
        
    for (int i = 0; i < numStates; i++)
        for (int j = 0; j < numStates; j++)
            q[i][j] *= scaleFactor;
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double sharingRates = 0.0;
    
    // CharMatrix(Tree* t, double** q, int ns, double* freqs, int nc, RandomVariable* rng, double alpha, double beta, std::vector<double> sharingTimes)
    charMatrix = new CharMatrix(tree, q, 2, f, 100, &rng, 9.0, 1.0, sharingRates);
    //charMatrix->print();
    [treeview setNeedsDisplay:YES];
}

// default values
- (id)init {
    if (self = [super init]) {
        [self setSpeciationRate:1.0];
        [self setExtinctionRate:0.1];
        [self setDuration:1.0];
    }
    return self;
}

- (Tree*)getTree {
    return tree;
}

@end
