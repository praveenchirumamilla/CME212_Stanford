/* Test script for the SDLViewer and Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point list
 *
 * Prints
 * A B
 * where A = number of nodes
 *       B = number of edges
 * and launches an SDLViewer to visualize the system
 */

#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  Graph<double> graph;
  std::vector<typename Graph<double>::node_type> nodes;

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

  // Print number of nodes and edges
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;


  // Launch a viewer
  CME212::SDLViewer viewer;
  viewer.launch();

  // add nodes using iterator
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map); 
  viewer.center_view(); 
  
  // Set the viewer
  //viewer.draw_graph_nodes(graph);
  //viewer.draw_graph(graph);
  //viewer.center_view();

  return 0;
}
