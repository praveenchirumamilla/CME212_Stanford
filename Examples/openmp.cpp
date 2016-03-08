#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>

#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP

#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/counting_iterator.h>
using thrust::make_permutation_iterator;
using thrust::make_counting_iterator;

#include <thrust/transform.h>
#include <thrust/system/omp/execution_policy.h>


#include "CME212/Util.hpp"

// RAII helper class for timing code
struct Timer {
  std::string msg;
  CME212::Clock clock;
  Timer(const std::string& s) : msg(s) {
    clock.start();
  }
  ~Timer() {
    double elapsed = clock.seconds();
    std::cout << msg << ": " << elapsed << "s" << std::endl;
  }
};


// Generic parallel-transform implementation that checks the iterator category
template <class Iter, class Function>
void parallel_transform(Iter first, Iter last, Function f) {
  // Check (at compile-time) the category of the iterator
  using category = typename std::iterator_traits<Iter>::iterator_category;
  static_assert(std::is_same<category, std::random_access_iterator_tag>::value,
                "parallel_transform requires random-access iterators!");

#pragma omp parallel for
  for (auto it = first; it < last; ++it)
    *it = f(*it);
}


int main()
{
  unsigned N = 100000000;
  std::vector<double> a;

  // Normal loop
  a = std::vector<double>(N,4);
  { Timer timer("Serial");
    for (auto it = a.begin(); it < a.end(); ++it)
      *it = std::exp(std::sqrt(*it));
  }

  // Parallel loop
  a = std::vector<double>(N,4);
  { Timer timer("OMP parallel for");
#pragma omp parallel for
    for (auto it = a.begin(); it < a.end(); ++it)
      *it = std::exp(std::sqrt(*it));
  }

  // Wrapped in a function
  a = std::vector<double>(N,4);
  { Timer timer("parallel transform");
    // Create a lambda function and call parallel_transform
    auto func = [](const double& ai) { return std::exp(std::sqrt(ai)); };
    parallel_transform(a.begin(), a.end(), func);
  }

  // Create a thrust iterator and call the function
  a = std::vector<double>(N,4);
  { Timer timer("parallel_transform Permutation Iter");
    // Make a trivial, but still random-access iterator pair
    auto begin = make_permutation_iterator(a.begin(),
                                           make_counting_iterator<int>(0));
    auto end   = make_permutation_iterator(a.end(),
                                           make_counting_iterator<int>(N));
    // Call the parallel_transform, still works even with weird iterators!
    auto func = [](const double& ai) { return std::exp(std::sqrt(ai)); };
    parallel_transform(begin, end, func);
  }

  /*
  // This won't compile!
  a = std::vector<double>(N,4);
  { Timer timer("Won't compile!");
    // Filter iterators are (at best) bi-directional, not random-access!
    auto f_true = [](double) { return true; };
    auto begin = boost::make_filter_iterator(f_true, a.begin(), a.end());
    auto end   = boost::make_filter_iterator(f_true, a.begin(), a.end());
    // Create a lambda function and call parallel_transform
    auto func = [](double ai) { return std::exp(std::sqrt(ai)); };
    parallel_transform(begin, end, func);
  }
  */

  // Use a thrust algorithm
  a = std::vector<double>(N,4);
  { Timer timer("thrust::transform");
    auto func = [](const double& ai) { return std::exp(std::sqrt(ai)); };
    thrust::transform(thrust::omp::par, a.begin(), a.end(), a.begin(), func);
  }

  return 0;
}
