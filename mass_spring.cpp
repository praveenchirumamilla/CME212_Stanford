/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"



// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
};

struct EdgeData {
  double L;
  double K;
};

// HW2 #1 YOUR CODE HERE
// Define your Graph type
typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;



/* sphere constraint */
struct sphereConstraint{

   Point c = Point(0.5, 0.5, -0.5);
   double r = 0.15;
   void operator()(GraphType& g, double){
    for(auto it = g.node_begin(); it != g.node_end(); ++it)
    {
	Node n = (*it);
	if(norm(n.position() - c) < r){
	    Point R = (n.position() - c)/norm(n.position() - c);
	    n.position() = c + r*R;
	    n.value().vel = n.value().vel - dot(n.value().vel, R)*R;
	}
     }
    }   
};

/* plane constraint */
struct planeConstraint{

  void operator()(GraphType& g, double){
    for(auto it = g.node_begin(); it != g.node_end(); ++it)
    {
	Node n = (*it);
	if(dot(n.position(), Point(0, 0, 1)) < -0.75){
	  n.position().elem[2] = -0.75;
	  n.value().vel.elem[2] = 0;
  	}
     }
   }
};

/* remove sphere constraint */

struct removeSphereConstraint{

   Point c = Point(0.5, 0.5, -0.5);
   double r = 0.15;
   void operator()(GraphType& g, double){
   auto it = g.node_begin();
   while(it != g.node_end()){
	if(norm((*it).position() - c) < r){
	  it = g.remove_node(it);
	}else{
	  ++it;
	}
    } 
   }  
};  

template<typename Con1, typename Con2>
struct combineTwoConstraints
{
 	Con1 constraint1;
	Con2 constraint2;
	
	combineTwoConstraints(Con1 c1 = Con1(), Con2 c2 = Con2()):constraint1(c1), constraint2(c2){}
	void operator()(GraphType& g, double){
		constraint1(g, 0);
		constraint2(g, 0);
	}
};

template<typename C1, typename C2>
combineTwoConstraints<C1, C2> totalConstraint(C1 c1 = C1(), C2 c2 = C2()){
	return combineTwoConstraints<C1, C2>(c1, c2);
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {

   /* I am getting an error with planeconstraint which i cannot able to resolve at this time */
  //auto tempConstraint = totalConstraint(removeSphereConstraint(), PlaneConstraint());
  
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    if(n.position() != Point(0,0,0) && n.position() != Point(1, 0, 0)){

      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }

    if(n.position() == Point(0,0,0) || n.position() == Point(1, 0, 0)){

      // set its velocity as 0
      n.value().vel = Point(0, 0, 0);
    }
  }
  

  // Compute the t+dt velocity
  //tempConstraint(g, 0);
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE

    if(n.position() == Point(0,0 ,0) || n.position() == Point(1, 0, 0)){
	return Point(0, 0, 0);
    }

    /* Spring Force & gravity force */
    Point spring = Point(0, 0, 0);
    Point gravity;
    gravity = Point(0, 0, -grav)*n.value().mass;

    /* sum up all spring forces */
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
	Edge e = *it;
	spring += (e.value().K*(e.node2().position() - e.node1().position())/e.length())*(e.length()- e.value().L); 
        //spring += (100*(e.node2().position() - e.node1().position())/e.length()*(e.length()- .1)); 
    }

    (void) t;     // silence compiler warnings
    return (spring + gravity);
  }
};

/* Gravity Force */
struct GravityForce{
   template <typename NODE>
   Point operator()(NODE n, double t){
	(void) t;
	return (Point(0, 0, -grav)*n.value().mass);
   }
};

/* Spring Force */	
struct MassSpringForce{

   template <typename NODE>
   Point operator()(NODE n, double t){
     Point spring = Point(0, 0, 0);
 
     /* sum up all spring forces */
     for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
	Edge e = *it;
        spring += (e.value().K*(e.node2().position() - e.node1().position())/e.length())*(e.length()- e.value().L);
     }

     (void) t;
     return spring;
   }
};

/* Damping force */
struct DampingForce	
{
   template <typename NODE>
   Point operator()(NODE n, double t){
	return (-(c*n.value().vel));
   }
   static double c;
};
double DampingForce::c = 0;

/* Combined force of F1 & F2 */
template<typename F1, typename F2>
struct CombinedForce{
	F1 force1;
	F2 force2;
        CombinedForce(F1 fo1=F1(), F2 fo2 = F2()):force1(fo1),force2(fo2){}

	template <typename NODE>
   	Point operator()(NODE n, double t)
        {
	  (void) t;
	  return (force1(n, 0) + force2(n, 0));
        }
};

/* combine two forces */
template<typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 fo1 = F1(), F2 fo2 = F2()){
  return CombinedForce<F1, F2>(fo1, fo2);
}

/* Combine three forces */
template<typename F1, typename F2, typename F3>
CombinedForce<F1, CombinedForce<F2, F3> > make_combined_force(F1 fo1, F2 fo2, F3 fo3){
	return make_combined_force(fo1, make_combined_force(fo2, fo3));
}


int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
//#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  /* Set mass and velocity */
  for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
	(*it).value().mass = float(1)/graph.size();
	(*it).value().vel = Point(0,0,0);
  }

  /* Set K & L */
  for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
     for(auto ij = (*it).edge_begin(); ij != (*it).edge_end(); ++ij){
	(*ij).value().K = 100;
	(*ij).value().L = (*ij).length();
     }
  } 

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0;
  double t_end = 5.0;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cout << "t = " << t << std::endl;
    //symp_euler_step(graph, t, dt, Problem1Force());
    symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()));

    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    /* for sphere removal constraint */
    viewer.clear();
    node_map.clear();

    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CME212::sleep(0.001);
  }

  return 0;
}
