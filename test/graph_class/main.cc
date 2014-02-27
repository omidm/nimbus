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
  * This file is written to check the functionality and performance of the
  * graph class.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <mcheck.h>
#include <stdlib.h>
#include <iostream> // NOLINT
#include "shared/nimbus.h"
#include "shared/nimbus_types.h"
#include "shared/graph.h"
#include "./vertex_entry.h"

#define VERTEX_NUM 5
#define GRAPH_NUM 1

int main(int argc, char *argv[]) {
  mtrace();

  nimbus::nimbus_initialize();

  int input;
  std::cout << "Start... ";
  std::cin >> input;

  for (int i = 0; i < GRAPH_NUM; i++) {
    nimbus::Graph<VertexEntry, int>* graph = new nimbus::Graph<VertexEntry, int>();
    std::cout << graph << std::endl;
    for (int j = 0; j < VERTEX_NUM; j++) {
      VertexEntry entry;
      graph->AddVertex(j, &entry);
    }
    for (int j = 0; j < VERTEX_NUM; j++) {
      for (int m = 0; m < j; m++) {
        graph->AddEdge(j, m);
      }
    }
    for (int j = 0; j < VERTEX_NUM; j++) {
      for (int m = 0; m < j; m++) {
        graph->AddEdge(j, m);
      }
    }

    std::cout << "After insertion... ";
    std::cin >> input;

    // pick either of the following to check the remove vertex and remove edge.
    for (int j = 0; j < VERTEX_NUM; j++) {
      graph->RemoveVertex(j);
      std::cout << "removed vertex "<< j << std::endl;
    }

    delete graph;

    std::cout << "After deletion... ";
    std::cin >> input;
  }


  std::cout << "After for loop... ";
  std::cin >> input;

  return 0;
}

