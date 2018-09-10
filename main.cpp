// Basic libraries
#include <iostream>
#include <fstream>
#include <chrono>
#include <map>
#include <iterator>


// Boost libraries
#include <boost/math/distributions.hpp>


// Type definitions
#include "Typedefs.h"


// User included libraries
#include "RandomGraphs.h"



template<class T>
void printVector(std::vector<T> V)
{
  for (auto s = V.begin(); s!=V.end(); ++s)
    std::cout << *s << std::endl;
}





template<class T> void write_to_file(T G, std::string fName);
template<> void write_to_file(WeightedGraph G, std::string fName)
{
  // Write to file
  std::ofstream file(fName, std::ofstream::out);
  boost::graph_traits<WeightedGraph>::edge_iterator ei, ei_end;
  
  for (boost::tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei) {
    
    file << get(boost::vertex_index, G, source(*ei, G)) << "\t"
	 << get(boost::vertex_index, G, target(*ei, G)) << "\t"
	 << get(boost::edge_weight, G, *ei)
	 << std::endl;
  }
  
  file.close();
}








template<class T> void write_labels(T G, std::string fName);
template<> void write_labels(WeightedGraph G, std::string fName)
{
  // Write to file
  std::ofstream file(fName, std::ofstream::out);
  boost::graph_traits<WeightedGraph>::vertex_iterator vi, vi_end;
  
  for (boost::tie(vi, vi_end) = boost::vertices(G); vi != vi_end; ++vi) {
    
    file << get(boost::vertex_index, G, *vi) << "\t"
	 << G[*vi].label << std::endl;
  }
  
  file.close();
}







int main(int, char*[])
{
  
  // Graph parameters
  int NV = 50, k = 2, m = 3, nEdges = static_cast<int>(NV * 0.2);
  std::vector<int> cluster_labels = {1, -1};
  std::vector<double> w = {10.0, 12.0};
  double w_ICC = 3.0, sdev = 1.0;
  
  


  /************************** BA Graph *****************************/
  WeightedGraph G_BA(NV);
  std::vector<int> labels(NV, -666);
  
  BarabasiAlbert(G_BA, NV, m, nEdges, w[0], sdev);
  //BAClusterModel(G_BA, k, NV / 2, m, nEdges, w, w_ICC, sdev, cluster_labels);
  
  
  // Assign labels
  boost::graph_traits<WeightedGraph>::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = boost::vertices(G_BA); vi != vi_end; ++vi) {
    
    labels[get(boost::vertex_index, G_BA, *vi)] = G_BA[*vi].label;
  }
  
  write_to_file(G_BA, "BA_graph.txt");
  //write_labels(G_BA, "BA_labels.txt");
  
  

  
  /************************** SBM Graph ****************************/
  std::vector<double> p = {0.5, 1.0};
  boost::numeric::ublas::matrix<double> Q(k, k);
  double etha = 2.0;
  double c_in = 10;
  double c_out = c_in / etha;
  Q(0,0) = c_in / static_cast<double>(NV);
  Q(0,1) = c_out / static_cast<double>(NV);
  Q(1,0) = c_out / static_cast<double>(NV);
  Q(1,1) = c_in / static_cast<double>(NV);
  
  
  WeightedGraph G_SBM(NV);
  stochastic_block_model(G_SBM, cluster_labels, p, Q, w, w_ICC, sdev);
  
  
  // Assign labels
  for (boost::tie(vi, vi_end) = boost::vertices(G_SBM); vi != vi_end; ++vi) {
    
    labels[get(boost::vertex_index, G_SBM, *vi)] = G_SBM[*vi].label;
  }
  
  write_to_file(G_SBM, "SBM_graph.txt");
  write_labels(G_SBM, "SBM_labels.txt");
  
  return 0;
}

