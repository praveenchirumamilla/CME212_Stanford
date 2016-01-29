/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>
#include <queue>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"


/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };

   /* Returns true if node1 is closer to p_ than node2. */
   template <typename NODE>
   bool operator()(const NODE& node1, const NODE& node2) const {
    /*if(std::distance(p_, node1.position()) < std::distance(p_, node2.position())){
	return true;
    }else{
   	return false;
    } */

    double dist1 = 0, dist2 = 0;
    
    dist1 =((node1.position().x - p_.x) * (node1.position().x - p_.x) +
            (node1.position().y - p_.y) * (node1.position().y - p_.y) +
 	    (node1.position().z - p_.z) * (node1.position().z - p_.z));

    dist2 =((node1.position().x - p_.x) * (node1.position().x - p_.x) +
            (node1.position().y - p_.y) * (node1.position().y - p_.y) +
 	    (node1.position().z - p_.z) * (node1.position().z - p_.z));

    return dist1 < dist2;
  }
};


/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int>& g, const Point& point) {
  // HW1 #4: YOUR CODE HERE
  //(void) g, (void) point;

  /* find a nearest node to the point */
  //MyComparator mc(point);
  
  unsigned int maxDistance = 0;
  std::queue<Point::size_type> nodeQueue;
  std::vector<bool> visited(g.size(), false);

  /* Node iterator for nearest node to the point */
  Graph<int>::node_iterator it = std::min_element(g.node_begin(), g.node_end(), MyComparator(point));

  /* to start with set values of all nodes to -1:
       i). for root node set val to 0 & visited state. */
  for(Graph<int>::node_iterator temp = g.node_begin(); temp != g.node_end(); ++temp){
	(*temp).value() = -1;
     }
     (*it).value() = 0;
     visited[(*it).index()] = true;

  /* BFS Algortithm */
  nodeQueue.push((*it).index());
  while(!nodeQueue.empty()){
	Point::size_type id = nodeQueue.front();
        nodeQueue.pop();

        /* iterate edges from this origin node */
        for(Graph<int>::incident_iterator iit = g.node(id).edge_begin(); iit != g.node(id).edge_end(); ++iit){
	      Point::size_type nextId = (*iit).node2().index();
              if(!visited[nextId]){
			/* adjust the value of next node */
                        visited[nextId] = true;
			g.node(nextId).value() = g.node(nextId).value()+1;
	                maxDistance = g.node(nextId).value();

			nodeQueue.push(nextId);
               }
             }
          }
          return maxDistance;
}



int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  viewer.launch();

  // HW1 #4: YOUR CODE HERE
  // Use shortest_path_lengths to set the node values to the path lengths
  // Construct a Color functor and view with the SDLViewer
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();
  
  return 0;
}
