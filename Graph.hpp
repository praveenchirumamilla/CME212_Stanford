#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/iterator/counting_iterator.h"

using namespace std;
using namespace thrust;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 *
 * @tparam V Type of node values
 */
template <typename V, typename E = char>
class Graph {

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Type of node values. */
  typedef V node_value_type;

  /** Type of edge values. */
  typedef E edge_value_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

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

  /** Type of node uids
      Indexing type of nodes vector, unique for each node 
      Value type of idx2uid[size_type] */
  typedef unsigned uid_type;

 private:

  /** node_element is an object containing pertinent information about nodes,
   *   namely, a point, a node value, and an index.
   */
  struct node_element {
    Point point;
    node_value_type value;
    size_type idx;
  };

  /** adj_element is an object containing pertinent information about edges
   *   that are incident to a given node, namely a second node and an edge value.
   */
  struct adj_element {
    uid_type n2;
    edge_value_type value;
    bool operator==(const uid_type& n) const { return n2==n; }
  };

  // node vector of Points, with xyz position, node value, and idx
  // indexed by uid
  std::vector<node_element> nodes;
  // node vector of node uids
  // indexed by idx
  std::vector<uid_type> idx2uid;
  // adjacency vector of vectors,
  // where vector i contains adj_elements of adjacent node's uid 
  // and incident edge's value, for all adjacencies/incidencies of node i 
  std::vector<std::vector<adj_element>> edges;
  

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
  class Node : private totally_ordered<Node> {
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph->nodes[id].point;
    }

     /** Return this node's position modifiably. */
    Point& position() {
      return graph->nodes[id].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    uid_type index() const {
      return graph->nodes[id].idx;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const node_type& n) const { 
      return graph == n.graph and id == n.id; 
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const node_type& n) const {
      return (graph < n.graph) or (graph == n.graph and id < n.id);
    }

    /* Retrieves this node's value (able to be modified) */
    node_value_type& value() {
      return graph->nodes[id].value;
    }

    /* Retrieves this node's value (not able to be modified) */
    const node_value_type& value() const {
      return graph->nodes[id].value;
    }

    /* Returns this node's degree (number of nodes it is adjacent to) */
    size_type degree() const { 
      return graph->edges[id].size();
    }

    /* Returns iterator of type @a incident_iterator, 
     * with iterator at the start of this node's adjacency list 
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph, id, 0);
    }

    /* Returns iterator of type @a incident_iterator, with 
     * iterator at one past the end of this node's adjacency list 
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph, id, degree());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Allow Node to access parent Graph object
    Graph* graph; 
    // Node's id
    uid_type id;
    /** Private Constructor */
    Node(const Graph* graph_, uid_type id_)
        : graph(const_cast<Graph*>(graph_)), id(id_) {
    }
  };

  /** Return the number of valid nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return idx2uid.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * optional @param[in] nval The new node's value (defaults to empty)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  node_type add_node(const Point& position, 
                     const node_value_type& nval = node_value_type()) {
    node_element new_node;
    new_node.point = position;
    new_node.value = nval;
    new_node.idx =  size();
    idx2uid.push_back(nodes.size());
    nodes.push_back(new_node);
    // Add an empty vector to edges to store this new nodes adjacencies
    std::vector<adj_element> new_adj;
    edges.push_back(new_adj);
    return Node(this, idx2uid[size()-1]);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const node_type& n) const {
    return (this == n.graph and n.index() < num_nodes());
  }

  /** Return the node with idx @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    return Node(this, idx2uid[i]);
  }

  /** Remove node @a n and all incident edges
   *
   * If n is in our graph, return 1
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - n.degree()
   *
   * Else (if n is not in our graph), @posts above are NOT upheld
   * and return is 0;
   *
   * Complexity: no more than O(n.degree()*d) operations, 
   *   where d = max degree of vertices adjacent to n
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) 
      return 0;
    // erase all edges incident to this node by updating adjacency list
    for(size_type i = 0; i < edges[n.id].size(); ++i) {
      uid_type n2_id = edges[n.id][i].n2;
      for(size_type j = 0; j < edges[n2_id].size(); ++j) {
        if (edges[n2_id][j].n2 == n.id) {
          edges[n2_id][j] = edges[n2_id].back(); //swap
          edges[n2_id].pop_back(); //pop
          break;
        }
      }
    }
    edges[n.id].clear();

    size_type old_idx = nodes[n.id].idx;
    // swap
    idx2uid[old_idx] = idx2uid.back();
    // update the index of node to be swapped
    nodes[idx2uid[old_idx]].idx = old_idx;
    // pop
    idx2uid.pop_back();
    return 1;
  }

  template<typename P>
  void onBoundaryConf(P p){
    for(auto it = node_begin(); it != node_end(); ++it){
	if(p(*it) > -2){
		(*it).value().boundary = true;
	}else{
		(*it).value().boundary = false;
	}
     }
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
    }

    /** Return one node of this Edge */   
    node_type node1() const {
      return Node(graph, n1);
    }

    /** Return other node of this Edge */
    node_type node2() const {
      return Node(graph, n2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const edge_type& e) const {
      bool local_eq = e.nmin() == nmin() and e.nmax() == nmax();
      return (graph == e.graph and local_eq);
    }

    /** Test whether this edge is less than @a e in a global order. 
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const edge_type& e) const {
      bool first = nmin() < e.nmin();
      bool second = nmin() == e.nmin() and nmax() < e.nmax();
      // Within a graph, edge with smallest smaller node is <
      // If smaller nodes are =, compare second nodes
      bool local_order = first or second;
      // Allow for global ordering, too
      return (graph < e.graph) or (graph == e.graph and local_order);
    }

    /** Return edge length: cartesian distance between n1 and n2 */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return edge value (modifialbe) */
    edge_value_type& value(){
      auto& adj_e = graph->edges[nmin()];
      auto it = std::find(adj_e.begin(), adj_e.end(), nmax());
      assert(it != adj_e.end());
      return (*it).value;
    } 

    /** Return edge value (not modifialbe) */
    const edge_value_type& value() const{
      auto& adj_e = graph->edges[nmin()];
      auto it = std::find(adj_e.begin(), adj_e.end(), nmax());
      assert(it != adj_e.end());
      return (*it).value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Allow Edge to access parent Graph object
    Graph* graph; 
    // id of the two nodes that edge is incident to
    uid_type n1;
    uid_type n2; 
    /** Private Constructor */
    Edge(const Graph* graph_, size_type n1_, size_type n2_)
        : graph(const_cast<Graph*>(graph_)), n1(n1_), n2(n2_) {
    }
    /** Return smaller node id */
    uid_type nmin() const {
      return std::min(n1, n2);
    }
    /** Return larger node id */
    uid_type nmax() const {
      return std::max(n1, n2);
    }
  }; 

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes())
   */
  size_type num_edges() const {                        
    size_type sum_edg = 0;
    for (size_type i = 0; i < idx2uid.size(); ++i) {
      sum_edg += edges[idx2uid[i]].size();
    }
    return sum_edg/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_edges()), hopefully less
   */
  edge_type edge(size_type i) const __attribute__ ((deprecated)) {  
    edge_iterator e_it = edge_begin();
    for (size_type j = 0; j < i; ++j) {
      ++e_it;
    }
    return *e_it;                    
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const node_type& a, const node_type& b) const { 
    if (edges[b.id].size() < edges[a.id].size()) 
      return has_edge(b, a); // optimize to only check smaller adj list
    for (size_type i = 0; i < edges[a.id].size(); ++i) {
      if (edges[a.id][i].n2 == b.id) {
        return true;
      }  
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge (with modified value) if it already exists.
   *
   * @pre @a a and @a b are distinct valid nodes of this graph  
   *
   * @param[in] a,b Nodes that are the endpts of edge to be added
   * @param[in] nval Optional edge value (defaults to empty edge_value_type)
   *
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   *
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().  <<<-explain that old value gets overwritten
   *       Else,                        new num_edges() == old num_edges() + 1.
   * @post Edge(a,b).value() = nval (If edge already exists, value gets modified, too) 
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type add_edge(const node_type& a,
                     const node_type& b, 
                     const edge_value_type& nval = edge_value_type()) { 
    // Ensure we only add the value to the half of adj vector with first node smaller
    if(b < a) {  
      add_edge(b, a, nval);
    }
    if(!has_edge(a, b)) {
      adj_element adj_a; 
      adj_element adj_b; 
      adj_a.n2 = a.id;
      adj_b.n2 = b.id;
      edges[a.id].push_back(adj_b);
      edges[b.id].push_back(adj_a);
    }
    edge_type new_edge = Edge(this, a.id, b.id);
    new_edge.value() = nval;
    return new_edge;
  }

  /** Removes edge given two nodes
   * 
   * If edge is in our graph, return 1
   * @post new num_edges() == old num_edges() - 1
   *
   * Else (if edge is not in our graph), @post above is NOT upheld
   * and return is 0
   *
   * Complexity: at most O(n1.degree() + n2.degree()) 
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!(has_edge(n1, n2)))
      return 0;
    // Remove from n2's adjacency list
    for (size_type j = 0; j < edges[n2.id].size(); ++j) {
      if (edges[n2.id][j].n2 == n1.id) {
        edges[n2.id][j] = edges[n2.id].back();
        edges[n2.id].pop_back();
        break;
      }
    }
    // Remove from n1's adjacency list
    for (size_type j = 0; j < edges[n1.id].size(); ++j) {
      if (edges[n1.id][j].n2 == n2.id) {
        edges[n1.id][j] = edges[n1.id].back();
        edges[n1.id].pop_back();
        break;
      }
    }
    return 1;
  }

  /** Removes a given edge
   * 
   * If @a e is in our graph, return 1
   * @post new num_edges() == old num_edges() - 1
   *
   * Else (if @a e is not in our graph), @post above is NOT upheld
   * and return is 0
   *
   * Complexity: at most O(n1.degree() + n2.degree()) 
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Removes edge at location of a given an edge iterator
   * 
   * @return The next valid edge iterator (i.e., in lower triangular half of adjacency list).
   *         If @a e_it == old edge_end() - 1, then return is new edge_end() (not dereferenceable)
   *
   * @post If *(@a e_it) is in our graph, new num_edges() == old num_edges() - 1
   *       Else (if *(@a e_it) is not in our graph), new num_edges() == old num_edges()
   *
   * Complexity: at most O(n1.degree() + n2.degree()) 
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*(e_it));
    e_it.fix();
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    idx2uid.clear();
    edges.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  /* class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    // Element type. 
    typedef Node value_type;
    // Type of pointers to elements. 
    typedef Node* pointer;
    // Type of references to elements. 
    typedef Node& reference;
    // Iterator category. 
    typedef std::input_iterator_tag iterator_category;
    // Difference between iterators 
    typedef std::ptrdiff_t difference_type;

    // Construct an invalid NodeIterator. 
    NodeIterator() {
    }

    // Returns Node corresponding to NodeIterator's current position (node id) 
    Node operator*() const {
      return Node(graph, graph->idx2uid[node_idx]);
    } 

    // Returns NoteIterator with position (node id) incremented by one
    // *
    // * @pre @a node_id < graph->num_nodes()
    // 
    NodeIterator& operator++() {
      ++node_idx;
      return *this;
    }

    // Compares two NodeIterators, declaring equality if they have
    // *  the same graph pointer and same position (node id)
    //
    bool operator==(const NodeIterator& ni) const { 
      return graph == ni.graph and node_idx == ni.node_idx;
    }

   private:
    // Allow Graph to access NodeIterator's private member data and functions.
    friend class Graph;
    // Allow NodeIterator to access parent Graph object
    Graph* graph;
    // NodeIterator's 'position' is node_id,
    // necessarily in the range 0 <= @a node_id < graph->num_nodes()
    size_type node_idx;
    // Private Constructor 
    // *
    // * @pre @a node_id must be a valid node index, 
    // *       i.e. must be in the range 0 <= @a node_id < graph->num_nodes() 
     //
    NodeIterator(const Graph* graph_, size_type node_idx_)
        : graph(const_cast<Graph*>(graph_)), node_idx(node_idx_) {
    }
  };  */
  
  struct idx2node{
     friend class Graph;
   
     __host__ __device__
     Node operator() (size_type id) {return myGraph->node(id);}
     const graph_type* myGraph;

     public:
       idx2node(const graph_type* graph){
         myGraph = const_cast<graph_type*>(graph);
      }
  };
 
  struct NodeIterator : thrust::transform_iterator<idx2node, thrust::counting_iterator<size_type>, Node>{
      // Import super class's constructor
      using NodeIterator::transform_iterator::transform_iterator;
      
      // Custom constructor
      NodeIterator(const graph_type* myGraph, size_type id=0)
              : NodeIterator :: transform_iterator(thrust::counting_iterator<size_type>(id), idx2node(myGraph)){}
  }; 

  /** Returns iterator of type @a node_iterator, with
   *  iterator position at the start of graph's node list 
   */
  node_iterator node_begin() const {
     return node_iterator(this, 0);
  }

  /** Returns iterator of type @a node_iterator, with
   *  iterator position at one past the end of graph's node list 
   */
  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
  }

  /** Remove node at location of given node iterator 
   *
   * @return n_it node_iterator to same location in idx2uid (i.e., new n_it == old n_it)
   *            If (has_node(*(old n_it)) && old n_it == old node_end() - 1),
   *                new n_it == new node_end() (i.e., can't be dereferenced)
   *            Else if (has_node(*(old n_it))),
   *                *(new n_it) == *(++(old n_it))
   *            Else, *(new n_it) == *(old n_it)
   *
   * @post If *(@a n_it) is in our graph, new num_nodes() == old num_nodes() - 1
   *                                      new num_edges() == old num_edges() - n.degree()
   *       Else (if n is not in our graph), new num_nodes() == old num_nodes()
   *                                        new num_edges() == old num_edges() 
   *
   * Complexity: no more than O(n.degree()*d) operations, 
   *   where d = max degree of vertices adjacent to n
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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

    // Construct an invalid EdgeIterator. 
    EdgeIterator() {
    }

    /* Returns Edge corresponding to EdgeIterator's current 
        "position" within @a graph's adjacency list */
    Edge operator*() const {
      return Edge(graph, outer, graph->edges[outer][inner].n2);
    }
    
    /** Returns EdgeIterator with position incremented by one new
     *  edge, skipping over edges with node1 > node2
     *
     * @pre EdgeIterator's current position cannot be the last new edge
             i.e., if @a outer == graph->edges.size() - 1, then it must be 
     *       that for some i, graph->edges[@a outer][@a inner + i] < @a outer 
     *
     * @post After incrementation, graph->edges[@a outer][@a inner] < @a outer 
     */
    EdgeIterator& operator++() { 
      ++inner;
      fix();
      return *this;
    }

    /** Compares two EdgeIterators, declaring equality if they have
     *  the same graph pointer and same position in @a graph's adjacency list
     */
    bool operator==(const EdgeIterator& ei) const {
      return graph == ei.graph and outer == ei.outer and inner == ei.inner;
    }

   private:
    // Allow Graph to access EdgeIterator's private member data and functions.
    friend class Graph;
    // Allow EdgeIterator to access parent Graph object
    Graph* graph;
    // EdgeIterator's 'position' is graph->edges[@a outer][@a inner],
    uid_type outer;
    size_type inner;
    /** Private Constructor 
     * 
     *  Calls @a fix() to ensure EdgeIterator is at a viable position 
     */
    EdgeIterator(const Graph* graph_, uid_type outer_, size_type inner_)
        : graph(const_cast<Graph*>(graph_)), outer(outer_), inner(inner_) {
          fix();
    }

    /** Ensures EdgeIterator's position is an edge where node1 < node2, fixing
     *   positions that do not match this by incrementing to the next valid position.
     *
     * @pre EdgeIterator's current position cannot be the last new edge
             i.e., if @a outer == graph->edges.size() - 1, then it must be 
     *       that for some i, graph->edges[@a outer][@a inner + i] < @a outer 
     *
     * @post After incrementation, graph->edges[@a outer][@a inner] < @a outer 
     */
    void fix() {
      while (outer < graph->edges.size()) {
        while (inner < graph->edges[outer].size()) {
          while (graph->edges[outer][inner].n2 < outer) {  
            return;
          }
          ++inner;
        }
        ++outer;
        inner = 0;
      }
      return;
    }

  };

  /** Returns iterator of type @a edge_iterator, with
   *  iterator position at the start of graph's adjacency list 
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }

  /** Returns iterator of type @a edge_iterator, with
   *  iterator position at one past the end of graph's adjacency list 
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges.size()-1, edges[edges.size()-1].size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>{
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

    /* Returns Edge corresponding to IncidentIterator's current 
        "position" within @a graph's adjacency list */
    Edge operator*() const{
      return Edge(graph, root, graph->edges[root][idx].n2);
    }

    /** Increments the index corresponding to @a root's adjacencies
     * 
     * @pre @a idx < @a graph->edges[@a root].size() - 1;
     */
    IncidentIterator& operator++() {
      ++idx;
      return *this;
    }

    /** Compares two IncidentIterators, declaring equality if they have
     *  the same graph pointer and same position in @a graph's adjacency list
     */
    bool operator==(const IncidentIterator& ii) const {
      return graph == ii.graph and root == ii.root and idx == ii.idx;
    }

   private:
     // Allow Graph to access IncidentIterator's private member data and functions.
    friend class Graph;
    // Allow IncidentIterator to access parent Graph object
    Graph* graph;
    // IncidentIterator's 'position' is graph->edges[@a root][@a idx],
    uid_type root;
    size_type idx;

    /** Private Constructor
     * 
     * @pre @a root, @a idx must be a viable position:
     *       @a idx < @a graph->edges[@a root].size();
     *
     */
    IncidentIterator(const Graph* graph_, uid_type root_, size_type idx_)
        : graph(const_cast<Graph*>(graph_)), root(root_), idx(idx_) {
    }
  };

};

#endif // CME212_GRAPH_HPP
