/**
 * @file subgraph.cpp
 * Test script for viewing a subgraph from our Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <iterator>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator
    : private equality_comparable<filter_iterator<Pred,It>> {
 public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It> self_type;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
    // HW1 #4: YOUR CODE HERE
  }

  /** Returns original instance of a particular iterator type 
    */ 
  value_type operator*() const{
	return *it_;
  }
  
  /** Increments the current iterator to the next.
    */ 
  self_type& operator++() {	
     it_ = ++it_;
     while (!p_(*it_) && it_ != end_) it_ = ++it_;
     return *this;
  }

  /** Checks if two iterators for equality 
    */ 
  bool operator==(const self_type& other) const {
     return it_ == other.it_ && end_ == other.end_;	
  }

  /** Returns the starting value of iterator
    */ 
  const It begin() const {
     return it_;
  }

  /** Returns the ending value of iterator
    */
  const It end() const {
        return end_;
  }

 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.

/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
		return n.position().x < 0;
  }
};

/** predicate for selecting nodes in a cube of length .7 
  *
  * returns true if a node falls with in the cube of size 0.7 
  */ 
struct cubePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
		return (n.position().x < .35 && n.position().x > -.35
		     && n.position().y < .35 && n.position().y > -.35
		     && n.position().z < .35 && n.position().z > -.35);
  }
};


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
  // Use the filter_iterator to plot an induced subgraph.
  
  /* Uncomments the following 3 lines to test the given predicate for half graph */
  //SlicePredicate sp;
  //auto it_begin = make_filtered<SlicePredicate, Graph<int>::NodeIterator>(graph.node_begin(), graph.node_end(), sp);
  //auto it_end = make_filtered<SlicePredicate, Graph<int>::NodeIterator>(graph.node_end(), graph.node_end(), sp);
  
  /* Currently executing the predicate i implemented */
  cubePredicate cp;
  auto it_begin = make_filtered<cubePredicate, Graph<int>::NodeIterator>(graph.node_begin(), graph.node_end(), cp);
  auto it_end = make_filtered<cubePredicate, Graph<int>::NodeIterator>(graph.node_end(), graph.node_end(), cp);
  
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(it_begin, it_end, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  return 0;
}

//--COMMENT
//--Great job! Everything works as specified!
//--Try to fix your indentation. One easy way to do this with vim/vi is to go to visual mode and type gg=G :)
//--END

//--grade10