/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL
struct IdentityMatrix{

	/** Constructor for IdentityMatrix */
	IdentityMatrix(unsigned int row): rows(row){
	}

	/** Compute the product of a vector with this identity matrix
          */
	//template <typename Vector>
	//Vector operator*(const Vector& x) const{
	//	return x;
	//}
        
        /** Helper function to perform multiplication. Allows for delayed
   	 *   evaluation of results.
	 *  Assign::apply(a,b) resolves to an assignment operation such as
	 *    a += b, a -= b, or a = b.
	 *  @pre @size(v) == size(w) */
	template <typename VectorIn, typename VectorOut, typename Assign>
	void mult(const VectorIn& v, VectorOut& w, Assign) const{
		for(unsigned int i = 0; i < rows; ++i){
			Assign::apply(w[i], v[i]);
		}
	}

	/** Matvec forwards to MTL's lazy mat_cvec_multiplier operator */
	template <typename Vector>
	mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
	operator*(const Vector& v) const {
		return {*this, v};
	}

	unsigned int size() const{
		return rows*rows;
	}
		
	unsigned int num_rows() const{
		return rows;
	}

	private:
		//Empty
		unsigned int rows;
};

/** The number of elemenets in the matrix. */
inline std::size_t size(const IdentityMatrix& A){
	return A.size();
}

/** The number of rows in the matrix */
inline std::size_t num_rows(const IdentityMatrix& A){
	return A.num_rows();
}

/** The number of columns in the matrix */
inline std::size_t num_cols(const IdentityMatrix& A){
	return A.num_rows();
}

/** Traits that MTL uses to determine properties of our IdentityMatrix. */
namespace mtl{
namespace ashape{
	
/** Define IdentityMatrix to be a non-scalar type. */
template <>
struct ashape_aux<IdentityMatrix> {
	typedef nonscal type;
};
}  // end namespace ashape

/** IdentityMatrix implements the Collection concept
  * with value_type and size_type */
template <>
struct Collection<IdentityMatrix> {
	typedef double value_type;
	typedef unsigned size_type;
};
}


int main()
{
  // HW3: YOUR CODE HERE
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver

  /**  Initialize Identity matrix */
  unsigned int row = 40;
  IdentityMatrix I(row);
  
  /** Set b such that x == 1 is the solution, start with x = 0 */
  itl::pc::identity<IdentityMatrix> P(I);
  mtl::dense_vector<double> x(row, 1), b(row);
  b = I * x;
  x = 0;

  itl::cyclic_iteration<double> it(b, 100, 1e-11, 0.0, 5);
  cg(I, x, b, P, it);
  std::cout << x << std::endl;
  
  return 0;
}
