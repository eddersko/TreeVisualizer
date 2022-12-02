//
//  Tree.hpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#ifndef Tree_hpp
#define Tree_hpp
#include <vector>
#include <string>

class RandomVariable;
class Node;

class Tree {
    
public:
    Tree(void) = delete;
    Tree(std::string newickString);
    Tree(double lambda, double mu, double t);
    Tree(double lambda, double mu, double t, double sharingRate, double alpha, double beta, double** q, int numStates);
    void initializeDownPassSequence(void);
    void print(void);
    std::vector<Node*>& getDownPassSequence(void) {return downPassSequence;}
    std::vector<Node*> nodesAtTime(double time);
    void addSharingEvents(RandomVariable* rng, double rate, std::vector<Node*>& sourceNodes);
    Node* getRoot(void) { return root; }
    
private:
    std::vector<Node*> nodes;
    Node* root;
    Node* addNode(void);
    Node* addNode(int idx);
    void passDown(Node* p);
    std::vector<Node*> downPassSequence;
    std::vector<Node*> downPassSequenceSubtree;
    void initializeCoordinates(void);
    void initializeCoordinates(double t);    
    void listNode(Node* p, int indent);
    std::vector<std::string> parseNewickString(std::string ns);
    void reindex(void);
};

#endif /* Tree_hpp */
