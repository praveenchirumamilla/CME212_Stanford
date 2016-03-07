//--comment
//--Late submission received but late day already used for HW0. No credit.
//--END
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

#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include "math.h"
//#include "BoundingBox.hpp"

#include "Graph.hpp"

struct NodeInfo
{
	double val;
	bool boundary;
};

// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
typedef Graph<NodeInfo, double> GraphType;  //<  DUMMY Placeholder

class GraphSymmetricMatrix
{
	public:
		GraphSymmetricMatrix(GraphType* graph):myGraph(graph){
		}

		/** Helper function to perform multiplication. Allows for delayed
   	 	*   evaluation of results.
	 	*  Assign::apply(a,b) resolves to an assignment operation such as
	 	*    a += b, a -= b, or a = b.
	 	*  @pre @size(v) == myGraph->size() */
		template <typename VectorIn, typename VectorOut, typename Assign>
		void mult(const VectorIn& v, VectorOut& w, Assign) const{
		  for(auto it = myGraph->node_begin(); it != myGraph->node_end(); ++it){
		    unsigned int i = (*it).index();
		    int tot = 0;

	       	    //check if i is on boundary:
		    if((*it).value().boundary){
			tot = v[i];   // A[i, j] == 1 if i is on boundary.
		    }else{
		        for(auto ij = (*it).edge_begin(); ij != (*it).edge_end(); ++ij){	
			  if(!(*ij).node2().value().boundary){
                            tot += v[(*ij).node2().index()]; // if j is on boundary
			  }
			}
		    //}
		      tot -= (*it).degree() * v[i]; // L[i, j]	
		    }
	            Assign::apply(w[i], tot);
		 }
		}

		/** Matvec forwards to MTL's lazy mat_cvec_multiplier operator */
		template <typename Vector>
		mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
		operator*(const Vector& v) const {
			return mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>(*this, v); 

		}

	//helper function for number of rows
	unsigned int rows() const{
		return myGraph->size();
	}

	private:
		GraphType* myGraph;

};

/** The number of elemenets in the matrix. */
inline std::size_t size(const GraphSymmetricMatrix& A){
	return A.rows()*A.rows();
}

/** The number of rows in the matrix */
inline std::size_t num_rows(const GraphSymmetricMatrix& A){
	return A.rows();
}

/** The number of columns in the matrix */
inline std::size_t num_cols(const GraphSymmetricMatrix& A){
	return A.rows();
}


/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  // HW3: YOUR CODE HERE
  auto it = g.node_begin();
  while(it != g.node_end()){
	if(bb.contains((*it).position())){
		g.remove_node(it);
	}else{
		++it;
	}
  }
  //(void) g; (void) bb;   //< Quiet compiler
  return;
}



int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
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

  /** Force Function */
  class forceFun{
    public:
	double operator()(GraphType::node_type n){
	  return 5*cos(norm_1(n.position()));
	}
  };

  /** Boundary Condition */
  class boundaryCondition{
  public:
	double operator()(GraphType::node_type n){
	  if(norm_inf(n.position()) == 1)
		return 0;
	  else if((norm_inf(n.position() - Point(0.6, 0.6, 0)) < 0.2) || (norm_inf(n.position() - Point(-0.6, -0.6, 0)) < 0.2))
		return -0.2;
          else if((norm_inf(n.position() - Point(-0.6, 0.6, 0)) < 0.2) || (norm_inf(n.position() - Point(0.6, -0.6, 0)) < 0.2))
		return -0.2;
	  else if(BoundingBox(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1)).contains(n.position()))		
		return 1;
	  else
		return -2;
	}
  };

  /** b Definition */
  forceFun f;
  boundaryCondition g;
  mtl::dense_vector<double> b(graph.size());

  for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
     if((*it).value().boundary){
	b[(*it).index()] = g(*it);
     }else{
	b[(*it).index()] = h*h*f(*it);
        for(auto ij = (*it).edge_begin(); ij != (*it).edge_end(); ++ij){
	  if((*ij).node2().value().boundary)
		b[(*it).index()] -= g((*ij).node2());
	}
     }
  }

  /** GraphSymmetricMatrix A using graph */
  GraphSymmetricMatrix A(&graph);

  /** SDLViewer */
  CME212::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map); 
  viewer.center_view();

  /** Solve Au =b */
  itl::pc::identity<GraphSymmetricMatrix> P(A);
  mtl::dense_vector<double> x(graph.size(), 0.0);
  
  itl::cyclic_iteration<double> iter(b, 100, 1e-11, 0.0, 10);
  cg(A, x, b, p, iter);

  return 0;
}

//--grade0
