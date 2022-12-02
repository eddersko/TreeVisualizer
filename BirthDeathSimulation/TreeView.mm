//
//  TreeView.m
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#import "TreeView.h"
#import "AppController.h"
#include "Tree.hpp"
#include "Node.hpp"
#include <vector>
#include <set>

@implementation TreeView

- (void)drawRect:(NSRect)dirtyRect {
    [super drawRect:dirtyRect];
    
    // color the background
    NSRect bounds = [self bounds];
    [[NSColor whiteColor] set];
    NSRectFill(bounds);
    
    Tree* t = [appController getTree];
    if (t != NULL) {
        std::vector<Node*>& dpseq = t->getDownPassSequence();
        for (Node* p : dpseq) {
            
            Node* pAncs = p->getAncestor();
            
            if (pAncs != NULL) {
                
                double x = p->getCoordinates().x * bounds.size.width;
                double y1 = p->getCoordinates().y * bounds.size.height;
                double y2 = pAncs->getCoordinates().y * bounds.size.height;
                
                NSPoint p1 = NSMakePoint(x, y1);
                NSPoint p2 = NSMakePoint(x, y2);
                [[NSColor blackColor] set];
                NSBezierPath* line = [NSBezierPath bezierPath];
                [line moveToPoint:p1];
                [line lineToPoint:p2];
                [line setLineWidth:1.0];
                [line stroke];
                
                
                // if source or dest not null
                
                if (p->getDest() != NULL) {
                    double x2 = p->getDest()->getCoordinates().x * bounds.size.width;
                    
                    NSPoint p1 = NSMakePoint(x, y1);
                    NSPoint p2 = NSMakePoint(x2, y1);
                    [[NSColor redColor] set];
                    NSBezierPath* line = [NSBezierPath bezierPath];
                    [line moveToPoint:p1];
                    [line lineToPoint:p2];
                    [line setLineWidth:1.0];
                    [line stroke];
                    
                    // draw arrowhead
                    
                    [line moveToPoint:p2];
                    
                    NSPoint p3;
                    NSPoint p4;
                    
                    if (x2 > x) {
                        p3 = NSMakePoint(x2 - (0.005*bounds.size.width), y1 - (0.005*bounds.size.height));
                        p4 = NSMakePoint(x2 - (0.005*bounds.size.width), y1 + (0.005*bounds.size.height));
                    } else {
                        p3 = NSMakePoint(x2 + (0.005*bounds.size.width), y1 - (0.005*bounds.size.height));
                        p4 = NSMakePoint(x2 + (0.005*bounds.size.width), y1 + (0.005*bounds.size.height));
                    }
                    
                    [line moveToPoint:p2];
                    [line lineToPoint:p3];
                    [line stroke];
                    [line moveToPoint:p2];
                    [line lineToPoint:p4];
                    [line stroke];
                }
                
               /* if (p->getSource() != NULL) {
                    double x2 = p->getSource()->getCoordinates().x * bounds.size.width;
                    
                    NSPoint p1 = NSMakePoint(x, y1);
                    NSPoint p2 = NSMakePoint(x2, y1);
                    [[NSColor redColor] set];
                    NSBezierPath* line = [NSBezierPath bezierPath];
                    [line moveToPoint:p1];
                    [line lineToPoint:p2];
                    [line setLineWidth:1.0];
                    [line stroke];
                }*/
                
            }
            
            std::set<Node*> pDescs = p->getDescendants();
            
            if (pDescs.size() != 0) {
                
                double y = p->getCoordinates().y * bounds.size.height;
                double x1 = p->getCoordinates().x * bounds.size.width;
                
                for (Node* n : pDescs) {
                    double x2 = n->getCoordinates().x * bounds.size.width;
                    NSPoint p1 = NSMakePoint(x1, y);
                    NSPoint p2 = NSMakePoint(x2, y);
                    [[NSColor blackColor] set];
                    NSBezierPath* line = [NSBezierPath bezierPath];
                    [line moveToPoint:p1];
                    [line lineToPoint:p2];
                    [line setLineWidth:1.0];
                    [line stroke];
                }
                
            }
        }
    }
}

@end
