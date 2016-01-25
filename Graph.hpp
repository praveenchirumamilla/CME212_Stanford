#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <vector>
#include <map>
#include <array>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

   /* Graph containers */
   vector<Point> nodes;
   //vector<std::array<int,2> > edges;
   vector<std::array<unsigned int,2> > edges;
 
public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
 	/* Nothing to define here */   
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
     /* Node's constructor is defined in Node's private section */ 
    }

    /** Return this node's position. */
    const Point& position() const {
	/* return the point */
	return myGraph->nodes[id]; 	 
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
    
      return this->id;
   }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      
      return (this->id == n.id) && (this->myGraph == n.myGraph); 
    }


    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    //--comment
    //--This does not preserve trichotomy
    //--START
    bool operator<(const Node& n) const {
    
      return (this->id < n.id) && (this->myGraph < n.myGraph); 
     }
    //--END
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	
    size_type id;
    graph_type *myGraph;

    Node(size_type id, const graph_type* myGraph){
	this->id = id;
	this->myGraph = const_cast<graph_type*> (myGraph);
    };

   };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    
    return nodes.size(); 
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
	
      nodes.push_back(position);
      Node n(size()-1, this); 
      return n; 

   }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
  
    return n.id < size();
   }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    /* return the node */
    Node n(i, this); 
    return n;
 }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      
      Node n(first_node, myGraph); 
      return n;      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      Node n(second_node, myGraph); 
      return n;      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      
	return ((((e.first_node == this->first_node) && (e.second_node == this->second_node)) || 
                ((e.first_node == this->second_node) && (e.second_node == this->first_node))) &&
		(e.myGraph == this->myGraph)); 

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    //--comment
    //--This does not preserve trichotomy
    //--START
    bool operator<(const Edge& e) const {

        return ((this->myGraph < e.myGraph)|| ((this->myGraph == e.myGraph) && ((this->first_node < e.first_node 
                && this->second_node < e.second_node) ||
                (this->second_node < e.first_node && this->first_node < e.second_node)))); 
    
    }
    //--END
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
   
    size_type first_node;
    size_type second_node; 
    graph_type *myGraph;

    Edge(size_type first_node, size_type second_node, const graph_type* myGraph){
		this->first_node = first_node;
		this->second_node = second_node;	
		this->myGraph = const_cast<graph_type*> (myGraph);
	};

 };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
  
	return edges.size();
   }
  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    
       Edge e(edges[i][0], edges[i][1], this);
       return e;  
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    Edge e(a.id, b.id, this);
    for(size_type i = 0; i < edges.size(); i++){
	if( (e.first_node == edges[i][0] && e.second_node == edges[i][1]) ||
	    (e.first_node == edges[i][1] && e.second_node == edges[i][0])){
		return true;
	} 
	
    }
    return false;
   }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
  
   Edge e(a.id, b.id, this);
   std::array<size_type, 2> temp = {a.id, b.id};  
   if(!has_edge(a, b)){
	edges.push_back(temp);	
    }
    return e;	
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
	
    nodes.clear();
    edges.clear();
  }

 private:

  // this is space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
//--comment
//--You can improve your choice of data structures
//--For example, vectors and pairs are a good idea,
//-- but don't give you O(1) addition of edges
//--END
//--grade9