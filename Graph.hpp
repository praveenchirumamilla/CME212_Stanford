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

template <typename V>
class Graph {
 private:

   /* Graph containers:
	i). I am using two individual vectors for storing nodes & values.
       ii). I implemented edges as both adjacency as well as vector of pairs.
      iii). I used Edge implementation with vector of pairs for implementing edge iterator. */
   vector<Point> nodes;
   vector< V > values;
   vector<std::array<unsigned int,2> > oldEdges;
   vector<vector<unsigned int> > edges;
 
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
  typedef V node_value_type; 

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

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
  class Node : private totally_ordered<Node>{
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
    bool operator<(const Node& n) const {
      return std::tie(this->myGraph, this->id) < std::tie(n.myGraph, n.id); 
     }

    /** Return this node's value. */
    node_value_type& value(){
   	return myGraph->values[id]; 
    }

    /** Return this node's value.
      *
      * Just a duplicate constant function
      */
    const node_value_type& value() const{
   	return myGraph->values[id]; 
    }

    /** Return this node's degree. */
    size_type degree() const{
	return myGraph->edges[this->id].size();
    }
   
   /** Return an incident_iterator for the first edge from this node. */
   incident_iterator edge_begin() const{
	//IncidentIterator i(this->id, edges[id].begin(), this->myGraph);
        IncidentIterator i(this->id, 0, this->myGraph);
	return i;
   }	

   /** Return an incident_iterator for the last edge from this node. */
   incident_iterator edge_end() const{
	//IncidentIterator i(this->id, edges[id].end(), this->myGraph);
	IncidentIterator i(this->id, myGraph->edges[id].size(), this->myGraph);
        return i;
   }	
 
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	
    size_type id;
    graph_type *myGraph;
	
    /** private constructor for node. */
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

   Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
      /* push the position and values into seperate vectors */ 
      nodes.push_back(position);
      values.push_back(val);
     
      /* push a temp vector when ever a node is created */	 
      vector<unsigned int> temp;     	
      edges.push_back(temp);  
      Node n(num_nodes()-1, this); 
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
  class Edge : private totally_ordered<Edge> {
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
    bool operator<(const Edge& e) const {
        return ((this->myGraph < e.myGraph)|| ((this->myGraph == e.myGraph) && ((this->first_node < e.first_node 
                && this->second_node < e.second_node) ||
                (this->second_node < e.first_node && this->first_node < e.second_node)))); 
    
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
   
    size_type first_node;
    size_type second_node; 
    graph_type *myGraph;

    /** private constructor for edge. */
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
	size_type tot_nodes = edges.size(); 
	size_type tot_edges = 0;
       
        for(unsigned int i = 0; i < tot_nodes; i++){
		tot_edges = tot_edges + edges[i].size();
	}
     
	/* I am counting every edge twice in the above logic,
	   so i half the count before i return */ 
	tot_edges /= 2; 
	return tot_edges;
   }

  /** Return the edge with index @a i.
    * @pre 0 <= @a i < num_edges()
    *
    * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
    */
  Edge edge(size_type i) const {
       size_type count = 0;
       size_type node1 = 0; 
       size_type node2 = 0;
    
       /* find the vector indices for ith edge */  
       for(size_type ni = 0; ni < edges.size(); ni++){
	for(size_type j = 0; j < edges[ni].size(); j++){
	 while(count <= i){
             if((edges[ni][j] <= ni)&& (j != edges[ni].size())){
		j++;
             } else if(j < edges[ni].size()){
		count += 1;
		j++;
	     }else if(j == edges[ni].size()){
		j = 0;
		ni++; 	
             }
	}
	node1 = ni;
        node2 = edges[ni][j-1];
	break;	
      }
      break; 
     }
     Edge e(node1, node2, this);
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
    size_type aDegree = edges[a.id].size();
    for(size_type i = 0; i < aDegree; i++){
	if(edges[a.id][i] == b.id){
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
	edges[a.id].push_back(b.id);	
	edges[b.id].push_back(a.id);

	/* Here i am also maintaining the vector of pairs for edge iterator. */
	oldEdges.push_back(temp);	
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

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    
   /** Returns a node instance for this index. */
   Node operator*() const{
        Node a(this->id, this->myGraph);
        return a;
    }
  
    /** Increments the iterator's index by one. */
    NodeIterator& operator++(){
	this->id = this->id + 1;
        return *this;
    } 
 
    /** Checks two iterators for equality. */ 
    bool operator == (const NodeIterator& x) const{
       return ((this->id == x.id) && (this->myGraph == x.myGraph));
    }
 
  private:
    friend class Graph;
   
    // HW1 #2: YOUR CODE HERE
    graph_type *myGraph;
    size_type id;
  
    /** private constructor for node iterator. */
    NodeIterator (size_type id, const graph_type* myGraph){
	this->id = id;
	this->myGraph = const_cast<graph_type*> (myGraph);
     };  
  };

  /** Returns the starting index for node iterator */
  node_iterator node_begin() const{
	NodeIterator a(0 , this); 
        return a;
  }

  /** Returns the end index for node iterator */
  node_iterator node_end() const{
        NodeIterator a(num_nodes() , this); 
	return a;  
  }
 
  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Returns an edge instance with this iterator's index */
    Edge operator*() const{
	Edge e(myGraph->oldEdges[this->id][0], myGraph->oldEdges[this->id][1], this->myGraph);
  	return e;	 
    }
    
    /** Increments the edge iterator's index by 1 */
    EdgeIterator& operator++(){
        this->id = this->id + 1;	
        return *this;
    }	

    /** Checks two edge iterators for equality. */
    bool operator == (const EdgeIterator& e) const{
       return ((this->id == e.id) && (this->myGraph == e.myGraph));
    } 

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type id; 
    graph_type *myGraph;
    
    /** private constructor for edge iterator. */
    EdgeIterator(size_type id, const graph_type* myGraph){
  		this->id = id;
		this->myGraph = const_cast<graph_type*> (myGraph);
 
    };
  };

  /** Returns the starting edge iterator instance */
  edge_iterator edge_begin() const{
       EdgeIterator e(0, this);
       return e;  
  }

  /** Returns the ending edge iterator instance */
  edge_iterator edge_end() const{
	size_type totEdges = oldEdges.size();
        EdgeIterator e(totEdges, this);
        return e;  
  }	 
 
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Returns instance of an edge with this iterator's index */
    Edge operator*() const{
	Edge e(this->firstNodeInd, myGraph->edges[this->firstNodeInd][this->secondNodeInd], this->myGraph); 
	return e;
    }

    /** Increments incident iterator's index by one */
    IncidentIterator& operator++(){
	 this->secondNodeInd = this->secondNodeInd + 1;	
         return *this;
    }

   /** Checks two incident iterators for equality */
   bool operator ==(const IncidentIterator& iit) const{
	return ((this->firstNodeInd == iit.firstNodeInd) && (this->secondNodeInd == iit.secondNodeInd) 
                 && (this->myGraph == iit.myGraph));
   } 
	
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type *myGraph;
    size_type firstNodeInd;
    size_type secondNodeInd;
  
    /** Private constructor for incident iterator */
    IncidentIterator (size_type firstNodeInd, size_type secondNodeInd, const graph_type* myGraph){
	this->firstNodeInd = firstNodeInd;
	this->secondNodeInd = secondNodeInd;
	this->myGraph = const_cast<graph_type*> (myGraph);
     };  
  };

 private:

  // this is space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
