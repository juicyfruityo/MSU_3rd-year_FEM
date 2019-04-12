#ifndef LOAD_MESH
#define LOAD_MESH

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>


// using node = 
struct node {
  int nid;
  float x, y;

  node(int nid, float x, float y)
      : nid(nid), x(x), y(y)  {}
};

std::vector<node> Nodes;

// using element = 
struct element {
  int eid;
  std::vector<int> num;
  std::vector<node> _node;

  element(int eid, int n1, int n2, int n3, int n4)
         :  eid(eid)
         {
           num.push_back(n1);
           num.push_back(n2);
           num.push_back(n3);
           num.push_back(n4);

           for (int i=0; i<4; ++i) {
               int j = 0;
               while (Nodes[j].nid != num[i])
                   ++j;
               _node.push_back(Nodes[j]);
               // можно потом совсем убрать Nodes из глобальной областиы
           }
         }
};

std::vector<element> Elements;


int load_mesh() {
  // std::vector<node> Nodes;
  // std::vector<element> Elements;
  std::ifstream f_nodes, f_elements;
  f_nodes.open("nodes_df_6x3.txt");
  f_elements.open("elements_df_6x3.txt");

  int nid, eid, n1, n2, n3, n4;
  float x, y;
  std::string line;
  std::size_t prev=0, next;

  while (std::getline(f_nodes, line)) {
      next = line.find(',', prev);
      nid = std::atoi(line.substr(prev, next-prev).c_str());
      prev = next + 1;
      next = line.find(',', prev);
      x = std::atof(line.substr(prev, next-prev).c_str());
      prev = next + 1;
      next = line.find(',', prev);
      y = std::atof(line.substr(prev, next-prev).c_str());
      // std::cout << nid << " " << x << " " << y << std::endl;
      prev = 0;

	  node tmp(nid, x, y);
	  Nodes.push_back(tmp);
      // Nodes.emplace_back(nid, x, y);
  }

  while (std::getline(f_elements, line)) {
      next = line.find(',', prev);
      eid = std::atoi(line.substr(prev, next-prev).c_str());
      prev = next + 1;
      next = line.find(',', prev);
      n1 = std::atoi(line.substr(prev, next-prev).c_str());
      prev = next + 1;
      next = line.find(',', prev);
      n2 = std::atoi(line.substr(prev, next-prev).c_str());
      prev = next + 1;
      next = line.find(',', prev);
      n3 = std::atoi(line.substr(prev, next-prev).c_str());
      prev = next + 1;
      next = line.find(',', prev);
      n4 = std::atoi(line.substr(prev, next-prev).c_str());
      prev = next + 1;
      next = line.find(',', prev);
      prev = 0;

	  element tmp(eid, n1, n2, n3, n4);
	  Elements.push_back(tmp);
      // Elements.emplace_back(eid, n1, n2, n3, n4);
  }

  f_nodes.close();
  f_elements.close();
  return 0;
}

#endif
