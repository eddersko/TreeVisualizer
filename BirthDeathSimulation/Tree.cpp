//
//  Tree.cpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#include "Tree.hpp"
#include "RandomVariable.hpp"
#include "Node.hpp"
#include <cmath>
#include <iostream>
#include "Msg.hpp"
#include "CognateSet.hpp"

/* Horizontal transfer events wih a simulated tree under a birth-death process
 
   lambda: birth rate
   mu: death rate
   t: duration/time
   sharingRate: rate of horizontal transfer
   alpha:
   beta:
   q: rate matrix
   numStates: number of states */
Tree::Tree(double lambda, double mu, double t, double sharingRate, double alpha, double beta, double** q, int numStates) {
    
    RandomVariable& rv = RandomVariable::randomVariableInstance();
    std::vector<Node*> activeLineages;
    std::vector<double> eq(2, 0.5);
    root = addNode();
        
    CognateSet* cs = new CognateSet(100, &rv, eq);
    root->setCognateSet(cs);
}

// Building tree object based on an existing tree topology
Tree::Tree(std::string newickString) {
    
    std::vector<std::string> tokens = parseNewickString(newickString);
    Node* p = NULL;
    bool readingBrlen = false;
    int idx = 0;
    
    for (int i = 0; i < tokens.size(); i++) {
        
        std::string token = tokens[i];
        
        if (token == "(") {
            
            Node* newNode = addNode();
            if (p == NULL) {
                root = newNode;
            } else {
                newNode->setAncestor(p);
                p->addDescendant(newNode);
            }
            p = newNode;
            
        } else if (token == ")" || token == ",") {
            
            if (p->getAncestor() == NULL)
                Msg::error("No ancestor found!");
            p = p->getAncestor();
            
        } else if (token == ":") {
            
            readingBrlen = true;
            
        } else if (token == ";") {
            
            if (p != root)
                Msg::error("Expecting p at root!");
            
        } else {
            if (readingBrlen == false) {
                
                Node* newNode = addNode();
                newNode->setAncestor(p);
                p->addDescendant(newNode);
                newNode->setName(token);
                newNode->setIndex(idx++);
                p = newNode;
                
            } else {
                double v = atof(token.c_str());
                p->setBrLen(v);
            }
            readingBrlen = false;
        }
    }
    initializeDownPassSequence();
    initializeCoordinates();
    
    for (Node* n : downPassSequence) {
        if (n->getNumDescendents() > 0)
            n->setIndex(idx++);
    }
    //print();
}

/* Simulating tree under birth-death process
    lambda: birth rate
    mu: death rate
    t: duration/time  */
Tree::Tree(double lambda, double mu, double t) {
    
    RandomVariable& rv = RandomVariable::randomVariableInstance();
    root = addNode(-1);
    Node* p = addNode(-1);
    root->setTime(0.0);
    root->addDescendant(p);
    p->setAncestor(root);
    std::vector<Node*> activeNodes;
    activeNodes.push_back(p);
    
    double currentTime = 0.0;
    while (currentTime < t && activeNodes.size() > 0) {
        
        double rate = activeNodes.size() * (lambda + mu);
        currentTime += -log(rv.uniformRv()) / rate;
        
        if (currentTime < t) {
            
            int whichNode = (int)(rv.uniformRv()*activeNodes.size());
            p = activeNodes[whichNode];
            p->setTime(currentTime);
            if (rv.uniformRv() < lambda / (lambda+mu)) { // speciation
                
                Node* q = addNode(-1);
                Node* r = addNode(-1);
                q->setAncestor(p);
                r->setAncestor(p);
                p->addDescendant(q);
                p->addDescendant(r);
                
                activeNodes[whichNode] = q; // replace p with tip node q
                activeNodes.push_back(r); // push r as current tip
                
            } else { // extinction
                
                // replace current node p with last node in array then pop last node from array
                activeNodes[whichNode] = activeNodes[activeNodes.size()-1];
                activeNodes.pop_back();
                
            }
        }
    }
    
    // set index and time for tip nodes
    for (int i = 0; i < activeNodes.size(); i++) {
        activeNodes[i]->setTime(t);
        activeNodes[i]->setIndex(i);
    }
    
    initializeDownPassSequence();
    initializeCoordinates(t);
    
    // set index for all internal nodes
    int idx = (int)activeNodes.size();    
    for (Node* p : downPassSequence) {
        if (p->getIndex()==-1)
            p->setIndex(idx++);
    }
    
    //print();
}

Node* Tree::addNode(void) {
    Node* newNode = new Node;
    nodes.push_back(newNode);
    return newNode;
}

Node* Tree::addNode(int idx) {
    Node* newNode = addNode();
    newNode->setIndex(idx);
    return newNode;
}

void Tree::initializeDownPassSequence(void) {
    downPassSequence.clear();
    passDown(root);
}

// tree traversal
void Tree::passDown(Node* p) {
    if (p != NULL) {
        std::set<Node*>& pDescendants = p->getDescendants();
        for (Node* n : pDescendants)
            passDown(n);
        downPassSequence.push_back(p);
    }
}

void Tree::initializeCoordinates(void) {
    
    double t = 0.0;
    for (Node* n : downPassSequence) {
        if (n->getNumDescendents()==0) {
            Node* p = n;
            double sum = 0.0;
            while (p->getAncestor() != NULL) {
                sum += p->getBrLen();
                p = p->getAncestor();
            }
            
            if (sum > t)
                t = sum;
        }
    }
    
    for (int i = (int)downPassSequence.size()-1; i >= 0; i--) {
        Node* p = downPassSequence[i];
        if (p->getAncestor()==NULL) {
            p->setTime(0.0);
        } else {
            p->setTime(p->getAncestor()->getTime()+p->getBrLen());
        }
    }
    
    int idx = 0;
    
    for (Node* n : downPassSequence) {
        if (n->getNumDescendents() == 0) {
            
            n->setXCoordinate(idx++);
            n->setYCoordinate(n->getTime()/t);
            
        } else {
            
            std::set<Node*>& nDescendants = n->getDescendants();
            double minX = 1000000.0;
            double maxX = 0.0;
            double maxY = 0.0;
            
            for (Node* m : nDescendants) {
                double x = m->getCoordinates().x;
                double y = m->getCoordinates().y;
                if (x < minX)
                    minX = x;
                if (x > maxX)
                    maxX = x;
                if (y > maxY)
                    maxY = y;
            }
            n->setXCoordinate((minX + maxX)*0.5);
            n->setYCoordinate(n->getTime()/t);
        }
    }

    idx--;
    for (Node* n : downPassSequence) {
        EKPoint c = n->getCoordinates();
        
        if (idx == 0) {
            n->setXCoordinate(0.5);
        } else {
            n->setXCoordinate(c.x/idx);
        }
    }
    
    
}

void Tree::initializeCoordinates(double t) {
    
    int idx = 0;
    
    for (Node* n : downPassSequence) {
        if (n->getNumDescendents() == 0) {
            
            n->setXCoordinate(idx++);
            n->setYCoordinate(n->getTime()/t);
            
        } else {
            
            std::set<Node*>& nDescendants = n->getDescendants();
            double minX = 1000000.0;
            double maxX = 0.0;
            double maxY = 0.0;
            
            for (Node* m : nDescendants) {
                double x = m->getCoordinates().x;
                double y = m->getCoordinates().y;
                if (x < minX)
                    minX = x;
                if (x > maxX)
                    maxX = x;
                if (y > maxY)
                    maxY = y;
            }
            n->setXCoordinate((minX + maxX)*0.5);
            n->setYCoordinate(n->getTime()/t);
        }
    }

    idx--;
    for (Node* n : downPassSequence) {
        EKPoint c = n->getCoordinates();
        
        if (idx == 0) {
            n->setXCoordinate(0.5);
        } else {
            n->setXCoordinate(c.x/idx);
        }
    }
}

void Tree::print(void) {
    listNode(root, 3);
}

void Tree::listNode(Node* p, int indent) {
    if (p != NULL) {
        std::set<Node*>& pDescendants = p->getDescendants();
        for (int i = 0; i < indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " ( ";
        for (Node* n : pDescendants)
            std::cout << n->getIndex() << " ";
        std::cout << ") ";
        if (p->getDest() != NULL)
            std::cout << "dest: " << p->getDest()->getIndex() << " ";
        if (p->getSource() != NULL)
            std::cout << "source: " << p->getSource()->getIndex() << " ";
        std::cout << p->getTime() << ", " << p->getBrLen() << " ";
        std::cout << "(" << p->getCoordinates().x << "," << p->getCoordinates().y << ") ";
        if (p == root)
            std::cout << "Root ";
        std::cout << std::endl;
        for (Node* n : pDescendants)
            listNode(n, indent+3);
    }
}


std::vector<std::string> Tree::parseNewickString(std::string ns) {
    
    std::vector<std::string> tokens;
    std::string word = "";
    
    for (int i = 0; i < ns.length(); i++) {
        char c = ns[i];
        
        if (c == '(' || c == ')' || c == ':' || c == ',' || c == ';') {
            
            std::string tok = std::string(1,c);
            if (word != "")
                tokens.push_back(word);
            word = "";
            tokens.push_back(tok);
            
        } else {
            
            std::string tok = std::string(1,c);
            word += tok;
            
        }
        
        
    }
    
#   if 0
    
    for (int i = 0; i < tokens.size(); i++)
        std::cout << i << "--" << tokens[i] << std::endl;
    
#   endif
    
    return tokens;

}

std::vector<Node*> Tree::nodesAtTime(double time) {
    
    std::vector<Node*> ndes;
    
    for (int i = (int)downPassSequence.size()-1; i >= 0; i--) {
        
        Node* p = downPassSequence[i];
        if (p->getAncestor() == NULL)
            p->setTime(0.0);
        else {
            //p->setTime(p->getAncestor()->getTime()+p->getBrLen());
            //std::cout << "p->getTime(): " << p->getTime() << ", p->getAncestor()->getTime(): " << p->getAncestor()->getTime() << std::endl;
            if ((p->getTime() >= time) && (p->getAncestor()->getTime() <= time)) {
                //std::cout << "time: " << time <<  ", p->getTime(): " << p->getTime() << ", p->getAncestor()->getTime(): " << p->getAncestor()->getTime() << std::endl;
                ndes.push_back(p);
            }
        }
        
    }
    
    return ndes;
}


void Tree::addSharingEvents(RandomVariable* rng, double rate, std::vector<Node*>& sourceNodes) {

    //print();
    // insert nodes into tree at source points
            
    for (Node* n : downPassSequence) {
        
        Node* nAncs = n->getAncestor();
        
        if (nAncs != NULL) {
            
            double v = nAncs->getTime();
            double len = n->getBrLen()+v;
            
            while (v < len) {
                
                nAncs = n->getAncestor();
                
                v += -log(rng->uniformRv())/rate;
                
                if (v < len) {
                    Node* p = addNode();
                    p->setTime(v);
                    //std::cout << p->getTime() << std::endl;
                    sourceNodes.push_back(p);
                                        
                    nAncs->addDescendant(p);
                    p->addDescendant(n);
                    nAncs->removeDescendant(n);
                    n->setAncestor(p);
                    p->setAncestor(nAncs);
                                            
                    Node* d = NULL;
                    p->setDest(d);
                }
            }
        }
    }
    
    initializeDownPassSequence();
    reindex();
                
    // add destination nodes to tree
    for (Node* n : sourceNodes) {
                
        double nTime = n->getTime();
        std::vector<Node*> activeNodes = nodesAtTime(nTime);
        //std::cout << "size of activeNodes: " << activeNodes.size() << std::endl;
        //if (activeNodes.size() == 1)         
          //  std::cout << "n's index is " << n->getIndex() << "and time is: " << nTime << std::endl;
        
        if (activeNodes.size() > 1) {
        
            Node* d = NULL;
            
            do {
                 d = activeNodes[(int)(rng->uniformRv()*activeNodes.size())];
            } while (d == n);
            
            Node* dAncs = d->getAncestor();
            
            if (dAncs == NULL)
                Msg::error("Why the hell is dAncs NULL?!");
            
            Node* dest = addNode();
            dest->setTime(nTime);
            
            d->setAncestor(dest);
            dest->setAncestor(dAncs);
            dAncs->removeDescendant(d);
            dAncs->addDescendant(dest);
            dest->addDescendant(d);
            
            dest->setSource(n);
            n->setDest(dest);
            //std::cout << "destination is: " << dest << std::endl;
            initializeDownPassSequence();
        }
    }
            
    // initialize branch lengths from time
    
    for (Node* n : downPassSequence) {
        
        Node* nAncs = n->getAncestor();
        if (nAncs != NULL)
            n->setBrLen(n->getTime()-nAncs->getTime());
        else
            n->setBrLen(0.0);
            
    }
    
    reindex();
    //print();
    
    
    
}

void Tree::reindex() {

    // check indices of tips
    std::set<int> tipIndices;
    int nT = 0;
        
    for (Node* n : downPassSequence) {
        if (n->getNumDescendents() == 0) {
            tipIndices.insert(n->getIndex());
            nT++;
        }
    }

    if (nT != tipIndices.size())
        Msg::error("Tip indices should be the same size as nT!");
    
    for (int i = 0; i < nT; i++) {
        std::set<int>::iterator it = tipIndices.find(i);
        if (it == tipIndices.end())
            Msg::error("Cannot find tip index!");
    }
    
    int interiorIdx = nT;
    for (Node* n : downPassSequence) {
        if (n->getNumDescendents() != 0) {
            n->setIndex(interiorIdx++);
        }
    }
}

/*
std::vector<Node*> desc;
desc.push_back(root);
int i = (int)downPassSequence.size()-1;
root->setIndex(i--);


while (desc.size() > 0) {
    Node* n = desc.back();
    desc.pop_back();
    std::set<Node*> descend = n->getDescendants();
    if (descend.size() != 0) {
        for (Node* m : descend) {
            m->setIndex(i--);
            desc.push_back(m);
        }
    }
}

*/
