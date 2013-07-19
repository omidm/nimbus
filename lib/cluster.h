#ifndef _CLUSTER
#define _CLUSTER

#include <iostream>
#include <string>
#include <set>

class Node;
class Link
{
  Node * ns;
  Node * ne;
  unsigned int capacity;
};

#define LinkSet std::set<Link>
#define LinkSetP std::set<Link*>

enum NodeType{ COMPUTER, SWITCH};

class Node
{
  public:
    std::string IP;
    int id;
    
    NodeType type;
    LinkSet linkSetP;
   
    Node();
    ~Node();
    virtual NodeType getType();
};

#define NodeSet std::set<Node>
#define NodeSetP std::set<Node *>

class Computer : public Node
{
  public:
    unsigned int memCap;
    unsigned int l1Cap;
    unsigned int l2Cap;
    unsigned int l3Cap;
    unsigned int coreNum;
    
    Computer();
    ~Computer();
    virtual NodeType getType();
    
};

class Switch : public Node
{
  public:
    unsigned int portNum;
    unsigned int crossSectBand;

    Switch();
    ~Switch();
    virtual NodeType getType();
};


class ClusterMap
{
  LinkSet linkSet;
  NodeSet nodeSet;

  void addNode(Node *);
  void delNode(Node *);

  void addLink(Link *);
  void delLink(Link *);

  unsigned int getLatency (Node *, Node *);
  unsigned int getCapacity (Node *, Node *);

  NodeSetP getRoute(Node *, Node *);

};





#endif
