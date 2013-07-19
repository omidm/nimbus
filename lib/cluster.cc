#include "cluster.h"

Node::Node()
{
};

Node::~Node()
{
};

NodeType Node::getType()
{
  return COMPUTER; 
};


Computer::Computer()
{
};

Computer::~Computer()
{
};

NodeType Computer::getType()
{
  return COMPUTER;
};

Switch::Switch()
{
};

Switch::~Switch()
{
};

NodeType Switch::getType()
{
  return COMPUTER;
};



