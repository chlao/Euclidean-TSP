#include "common.h"
#include "point.h"
#include "MST.h"
#include "Minmatching/PerfectMatching.h"

/*
This project is a starter code and wrappers for CSE101W15 Implementation project.

point.h - uniform random pointset generator

MST.h - minimum spanning tree

PerfectMatching.h - interface to min cost perfect matching code 

-------------------------------------
PerfectMatching is from the paper:

Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.

sourcecode : pub.ist.ac.at/~vnk/software/blossom5-v2.05.src.tar.gz

*/

void LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, float** adjacentMatrix, int N) {
	int e = 0;
	node_num = N;
	edge_num = N*(N-1)/2 ; //complete graph

	edges = new int[2*edge_num];
	weights = new int[edge_num];

	for(int i = 0; i < N ; ++i) {
		for(int j = i+1 ; j< N ; ++j) {
			edges[2*e] = i;
			edges[2*e+1] = j;
			weights[e] = adjacentMatrix[i][j];
			e++;
		}
	}

	if (e != edge_num) { 
		cout<<"the number of edge is wrong"<<endl;

		exit(1); 
	}
}

void PrintMatching(int node_num, PerfectMatching* pm) {
	int i, j;

	for (i=0; i<node_num; i++) {
		j = pm->GetMatch(i);
		if (i < j) printf("%d %d\n", i, j);
	}
}

int main() {
	set< pair<int,int> > generatedPointset;
	float** adjacentMatrix;
	int W, H, N;
	Point pointset;
	W = 11390;
	H = 11232;
	N = 8586;

	//N = 1000; 


	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl;

	pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation
	//pointset.printPointset();

	generatedPointset = pointset.getPointset();
	adjacentMatrix = pointset.getAdjacentMatrix();

	//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal
	MST mst(adjacentMatrix, N);
	mst.makeTree();
	//mst.printMST();

	//Deliverable B: Find TSP2 path from the constructed MST
	//You won't need any wrappers for B.

  	mst.makeTSP2(); 


	//Deliverable C: Find TSP1.5 path from the constructed MST
	//mst.makeTSP1_5(); 

	return 0;
}
