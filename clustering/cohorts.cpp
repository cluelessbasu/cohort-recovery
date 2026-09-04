// This code tests the decomposition function in the file Escape/Decomposition.h.
// 
// 
// C. XXXX-2, Aug 2023
// 

#include "Escape/GraphIO.h"
#include "Escape/Graph.h"
#include "Escape/Decomposition.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <inttypes.h>
#include <sstream>
#include <iomanip>


using namespace std;
using namespace Escape;

void printDecomposition(vector<Cluster> decomposition, CGraph* cg, string name)
{
  cout << "\n--------------------------------------------\n";
  cout << "Generating details of "+name+" decomposition\n";

  FILE* f = fopen((name+"-decomposition.txt").c_str(),"w");
  if (!f)
  {
      printf("could not write to output to decomp.txt\n");
      return;
  }

  FILE* fs = fopen((name+"-stats.txt").c_str(),"w");
  if (!fs)
  {
      printf("could not write to output to stats.txt\n");
      return;
  }

  // write out the decomposition and stats
  Cluster* clus_array = new Cluster[2*cg->nVertices+1];
  //cout<<"Here\n";
  //cout << decomposition.size();
  Count total_clustered = 0;
  Count* freq = new Count[cg->nVertices+1];
  for (VertexIdx i=0; i < cg->nVertices; i++)
      freq[i] = 0;
  Count maxsize = 0;

  // let's populate the array of clusters
  Count num_clusters = 0;
  //cout<<"Here\n";
  for (auto iter = begin(decomposition); iter != end(decomposition); iter++) // iterate over the decomposition
  {
      //cout<<"Here\n";
      clus_array[num_clusters] = *iter; // copy the cluster into the array
      num_clusters++; // increment the number of clusters
  }

  sort(clus_array, clus_array+num_clusters, clusterCompareDecreasing); // sort the clusters in decreasing order of size
  //cout<<"Here\n";

  Count total_internal = 0; // keep track of the cut edges
  //FILE* fcut = fopen((name+"-cut.txt").c_str(),"w");
  unordered_set<VertexIdx> total_vertices;
  for (Count i=0; i < num_clusters; i++)
  {
      Cluster next_clus = clus_array[i]; // get the next cluster
      unordered_set<VertexIdx> next = next_clus.vertices; // get the vertices of the next cluster
      for (auto iterc = begin(next); iterc != end(next); iterc++)
          fprintf(f,"%" PRId64 " ",*iterc); // printing the vertices of the cluster
      fprintf(f,"\n"); // putting a new line

      // now we write out the stats of the decomposition
      if (next_clus.nVertices > maxsize)
          maxsize = next_clus.nVertices; // keep track of the largest cluster seen thus far
      freq[next_clus.nVertices]++; // increment the number of clusters of the given size

      double density = 0; // variable for density
      if (next_clus.vertices.size() == 1) // singleton cluster
      {
          density = 0;
      }
      else
      {   
          density = 2*(next_clus.nEdges)/((double)(next_clus.nVertices)*(double)(next_clus.nVertices-1));
          total_vertices.insert(next_clus.vertices.begin(), next_clus.vertices.end());
          //total_clustered += next_clus.nVertices; 
          total_internal += next_clus.nEdges; 
          fprintf(fs,"%" PRId64 " %.2f\n",next_clus.nVertices,density);
      }
      // now, we print out the cut edges
      /*
      for (auto iterc = begin(next); iterc != end(next); iterc++)
      {
          VertexIdx vert = *iterc; // getting the next vertex of the cluster
          for (EdgeIdx j = cg->offsets[vert]; j < cg->offsets[vert+1]; j++) // loop over the neighbors of vert
          {
              VertexIdx nbr = cg->nbors[j]; // get the neighbor of vert
              if (nbr < vert) // tie breaking to ensure edge is printed only once
                  continue;
              if (next_clus.vertices.find(nbr) == next_clus.vertices.end()) // neighbor is not in cluster
                  fprintf(fcut,"%" PRId64 " %" PRId64 "\n",vert,nbr); // write down the edge
          }
      }*/

  }
  total_clustered = total_vertices.size();
  printf("Non-trivial clustered vertices: %" PRId64 "\n",total_clustered);
  cout << "\n--------------------------------------------\n";

  FILE* ffreq = fopen((name+"-freq.txt").c_str(),"w");
  for (Count i=0; i < maxsize+1; i++)
  {
      fprintf(ffreq,"Size = %" PRId64 ", freq = %" PRId64 ", total vertices = %" PRId64 "\n",i,freq[i],i*freq[i]);
  }

  delete[] clus_array; // free memory
  delete[] freq;

  fclose(fs);
  fclose(ffreq);
  //fclose(fcut);
  fclose(f);

}

void printScores(std::map<VertexIdx,double> scores, string name) 
{
  cout << "\n--------------------------------------------\n";
  printf("generating singleton scores\n");
  FILE* f = fopen((name+"-scores.txt").c_str(),"w");
  if (!f)
  {
      printf("could not write to scores.txt\n");
      return;
  }
  fprintf(f,"Vertex,Score\n");
  for (const auto& entry : scores) {
        fprintf(f, "%" PRId64 ",%f\n", entry.first, entry.second);
    }
  cout << "\n--------------------------------------------\n";
}





int main(int argc, char *argv[])
{
  Graph g;

  if(argc <= 3) {
    printf("Error: Missing arguments\n");
    exit(1);
  }
  string name(argv[2]);

  // get epsilon as an argument
  double eps = atof(argv[3]); // convert to double
  printf("eps = %f\n",eps);

  //get the mode as an argument
  char mode;
  if(argc <= 4) {
    printf("Error: Missing mode, defaulting to simple\n");
    mode = 's';
  } else {
    mode = argv[4][0];
  }
  printf("mode = %c\n",mode); //modes available: s simple, e expand

  int threshold = 0;
  if (mode == 'e') {
    // get the threshold as an argument
    if(argc <= 5) {
        printf("Error: Missing threshold value\n");
        exit(1);
    }
    threshold = atoi(argv[5]);
    printf("threshold = %i\n",threshold);
  }

  if (loadGraph(argv[1], g, 1, IOFormat::escape))
    exit(1);

  printf("Loaded graph\n");
  CGraph cg = makeCSR(g);

  // printing the graph
//   FILE* fg = fopen("graph.txt","w");
//   cg.print(fg);

  printf("Vertices = %" PRId64 "\n",cg.nVertices);
  printf("Edges = %" PRId64 "\n",cg.nEdges/2);
  
  std::stringstream ss;
  ss << std::fixed << std::setprecision(1) << eps; // converts 0.1 to "0.1"
  std::string eps_str = ss.str();

  // 2. Find the decimal point and replace it with 'p'
  size_t decimal_pos = eps_str.find('.');
  if (decimal_pos != std::string::npos) {
      eps_str.replace(decimal_pos, 1, "p");
  }

  vector<Cluster> decomposition;
  std::map<VertexIdx, VertexIdx> cluster_map;
  //decomposition = OverlappingExtract(&cg, eps, cluster_map); // get the decomposition
  //printDecomposition(decomposition,&cg,"cohorts-overlapping-seq-" + eps_str + "-" + name);
  decomposition = IndExtract(&cg, eps, cluster_map); // get the decomposition
  printDecomposition(decomposition,&cg,"cohorts-overlapping-mis-" + eps_str + "-" + name);


}
