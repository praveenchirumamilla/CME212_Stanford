/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"
#include "math.h"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


struct NodeData
{
  double val;   
  bool boundary;  
};

// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
typedef Graph<NodeData,double> GraphType;
class GraphSymmetricMatrix
{
  public:
  GraphSymmetricMatrix(GraphType* graph):myGraph(graph) {}

   /**  Helper function for multiplicayion 
     *  Assign::apply(p,q) resolves assignment operations like += -= =
     **/
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v, VectorOut& w, Assign) const{
    for (auto it = myGraph->node_begin(); it != myGraph->node_end(); ++it)
    {
      unsigned int i = (*it).index();
      double tot = 0;

      if ((*it).value().boundary){
        tot = v[i];
      }
      else{
        for (auto ij = (*it).edge_begin(); ij != (*it).edge_end(); ++ij){
          if (!(*ij).node2().value().boundary){
            tot += v[(*ij).node2().index()]; //  neiter i nor j is on boundary, and i!= j
          }
        }
        tot -= (*it).degree() * v[i];  // -deg(ni) when i is not on boundary
      }
      Assign::apply(w[i], tot);
    }
  }

  /** Matrix - vector multiplication */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector> operator*(const Vector& v) const {
    return mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>(*this, v);
  }

  //helper function for number of rows
  unsigned int rows() const{
    return myGraph->size();
  }

  private:
  GraphType* myGraph;
};

/** The number of elements in the matrix*/
inline std::size_t size(const GraphSymmetricMatrix& A){
  return A.rows()*A.rows();
}
/** The number of rows in the matrix*/
inline std::size_t num_rows(const GraphSymmetricMatrix& A){
  return A.rows();
}
/** The number of columns in the matrix*/
inline std::size_t num_cols(const GraphSymmetricMatrix& A){
  return A.rows();
}

/** Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl{
namespace ashape{
  /**Define GraphSymmetricMatrix to be a non-svalar type*/
template<>
struct ashape_aux<GraphSymmetricMatrix>
{
  typedef nonscal type;
};
}  //end namespace ashape

/**  GraphSymmetricMatrix implements the Collection concept
 *  with value_type and size_type*/
template<>
struct Collection<GraphSymmetricMatrix>
{
  typedef double value_type;
  typedef unsigned size_type;
};
}//end namesapce mtl


/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  // HW3: YOUR CODE HERE
  auto it = g.node_begin();
  while (it != g.node_end()){
    if (bb.contains((*it).position()))
    {
      g.remove_node(it);
    }
    else
      ++it;
  }

  //(void) g; (void) bb;   //< Quiet compiler
  return;
}

/** Functor to node position */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& node) {
    Point position = node.position();
    position.elem[2] = node.value().val;
    return position;
  }
};

/** Functor for node color */
struct NodeColor{
  CME212:: Color operator() (GraphType::node_type node){
    if (node.value().val < -1)
      return CME212::Color::make_heat(1);
    else if (node.value().val < 0.1)
      return CME212::Color::make_heat(0.2);
    else if (node.value().val < 1)
      return CME212::Color::make_heat(0.9);
    else
      return
        CME212::Color::make_heat(1);
  }
};


template <typename Real, typename PointFn, typename ColorFn>
class visual_iteration : public itl::cyclic_iteration<Real> 
{
   typedef itl::cyclic_iteration<Real> super;
   typedef visual_iteration self;
 public:

   template <class Vector>
   visual_iteration(GraphType* graph, CME212::SDLViewer* viewer, mtl::dense_vector<double>* x, const Vector& r0, int max_iter_, Real tol_, Real atol_ = Real(0), int cycle_ = 100)
     : super(r0, max_iter_, tol_, atol_, cycle_), myGraph(graph), viewer(viewer), x(x)
   {
   }
   

   bool finished() { return super::finished(); }

   template <typename T>
   bool finished(const T& r) 
   {
       //process the job and update each node's value.
       bool ret= super::finished(r);
       for (auto it = myGraph->node_begin(); it != myGraph->node_end(); ++it){
        (*it).value().val = (*x)[(*it).index()];
       }

       //update viewer
       auto node_map = viewer->empty_node_map(*myGraph);
       viewer->clear();
       viewer->add_nodes(myGraph->node_begin(), myGraph->node_end(), ColorFn(), PointFn(), node_map);
       viewer->add_edges(myGraph->edge_begin(), myGraph->edge_end(), node_map);
       return ret;
   }
 
   private:
   GraphType* myGraph;  
   CME212::SDLViewer* viewer;  
   mtl::dense_vector<double>* x;  
};


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);

  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }
 
  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.
  
  //boundary checks
  class Boundary{
  public:
    double operator() (GraphType::node_type n){
      if (norm_inf(n.position()) == 1)
        return 0;
      else if ((norm_inf(n.position()-Point(0.6, 0.6, 0)) < 0.2) || (norm_inf(n.position()-Point(-0.6, -0.6, 0)) < 0.2))
        return -0.2;
      else if ((norm_inf(n.position()-Point(-0.6, 0.6, 0)) < 0.2) || (norm_inf(n.position()-Point(0.6, -0.6, 0)) < 0.2))
        return -0.2;
      else if (Box3D(Point(0.6,0.2,1), Point(-0.6,-0.2,-1)).contains(n.position()))
        return 1;
      else
        return -5;  //jus non qualifier for a node which is not on the boundary
    }
  };

  //Forcing function in the Poisson equation
  class Force{
  public:
    double operator() (GraphType::node_type n){
      return 5*cos(norm_1(n.position()));
    }
  };

  // b, f , g functions
  Boundary g;
  Force f;
  mtl::dense_vector<double> b(graph.size());
  graph.onBoundaryConf<Boundary>(g);

  //configure  b
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    if ((*it).value().boundary){
      b[(*it).index()] = g(*it);   
    }else{
      b[(*it).index()] = h*h*f(*it);  
      for (auto ij = (*it).edge_begin(); ij != (*it).edge_end(); ++ij){
        if ((*ij).node2().value().boundary)
          b[(*it).index()] -= g((*ij).node2());
      }
    }
  }

  // Construct the GraphSymmetricMatrix A using the graph
  GraphSymmetricMatrix A(&graph);

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  // Solve Au = b using MTL.
  itl::pc::identity<GraphSymmetricMatrix> P(A);
  mtl::dense_vector<double> x(graph.size(), 0.0);
  visual_iteration<double, NodePosition, NodeColor> iter(&graph, &viewer, &x, b, 1000, 1.e-10, 0.0, 5);

  cg(A, x, b, P, iter);

  return 0;
}


