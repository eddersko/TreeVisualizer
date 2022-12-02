//
//  Node.hpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#ifndef Node_hpp
#define Node_hpp
#include <set>
#include <string>
#include <vector>
class CognateSet;
class RandomVariable;

struct EKPoint {
    double x;
    double y;
};

class Node {
  
    public:
        Node(void);
        void addDescendant(Node* p) { descendants.insert(p); }
        void removeDescendant(Node* p) { descendants.erase(p); }
        void removeAllDescendants(void) { descendants.clear(); }
        int getNumDescendents(void) { return (int)descendants.size(); }
        void setAncestor(Node* p) {ancestor = p; }
        Node* getAncestor(void) { return ancestor; }
        void setIndex(int i) {index = i;}
        int getIndex(void) {return index;}
        void setTime(double t) {time = t;}
        double getTime(void) {return time;}
        std::set<Node*>& getDescendants(void) {return descendants;}
        void setYCoordinate(double y) {coordinate.y = y;}
        void setXCoordinate(double x) {coordinate.x = x;}
        EKPoint getCoordinates(void) {return coordinate;}
        void setName(std::string n) {name = n;}
        std::string getName(void) {return name;}
        void setBrLen(double brLen) {v = brLen;}
        double getBrLen(void) {return v;}
        //double getTotalBrLen(void);
        int getState(void) {return state;}
        void setState(int s) {state = s;}
        void setCognateSet(CognateSet* cs) { myCognate = cs; }        
        CognateSet* getCognateSet(void) { return myCognate; }
        void setSource(Node* s) { source = s; }
        Node* getSource(void) { return source; }
        void setDest(Node* d) { destination = d; }
        Node* getDest(void) { return destination; }
    
    private:
        std::set<Node*> descendants;
        Node* ancestor;
        Node* destination;
        Node* source;
        std::string name;
        int index;
        double time;
        double v;
        EKPoint coordinate;
        int state;
        CognateSet* myCognate;
};

#endif /* Node_hpp */
