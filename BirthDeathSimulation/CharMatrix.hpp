//
//  CharMatrix.hpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/26/22.
//

#ifndef CharMatrix_hpp
#define CharMatrix_hpp
class Tree;
class RandomVariable;
class Node;

class CharMatrix {
    
public:
    CharMatrix(void) = delete;
    ~CharMatrix(void);
    CharMatrix(Tree* t, double** q, int ns, double* freqs, int nc, RandomVariable* rng);
    CharMatrix(double** q, int ns, double* freqs, int nc, RandomVariable* rng, double alpha, double beta, double sharingTimes);
    CharMatrix(Tree* t, double** q, int ns, double* freqs, int nc, RandomVariable* rng, double alpha, double beta, double sharingTimes);
    void print(void);
private:
    int numTaxa;
    int numChar;
    int numStates;
    int** matrix;
    double* resilience;
    void simulateSubTree(Node* n, Node* r, double** q, RandomVariable* rng, int site);
};

#endif /* CharMatrix_hpp */
