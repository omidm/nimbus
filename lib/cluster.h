/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

 /*
  * Nimbus abstraction of computational resources. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#ifndef NIMBUS_LIB_CLUSTER_H_
#define NIMBUS_LIB_CLUSTER_H_

#include <iostream>
#include <string>
#include <set>

class Computer;
typedef std::set<Computer*> Hosts;

class Node;
class Link {
  Node * ns;
  Node * ne;
  unsigned int capacity;
};

#define LinkSet std::set<Link>
#define LinkSetP std::set<Link*>

enum NodeType {COMPUTER, SWITCH};

class Node {
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

class Computer : public Node {
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

class Switch : public Node {
  public:
    unsigned int portNum;
    unsigned int crossSectBand;

    Switch();
    ~Switch();
    virtual NodeType getType();
};


class ClusterMap {
  LinkSet linkSet;
  NodeSet nodeSet;

  void addNode(Node * n);
  void delNode(Node * n);

  void addLink(Link * l);
  void delLink(Link * l);

  unsigned int getLatency(Node * n1, Node * n2);
  unsigned int getCapacity(Node * n1, Node * n2);

  NodeSetP getRoute(Node * n1, Node * n2);
};

#endif  // NIMBUS_LIB_CLUSTER_H_


