#include "common.h"
#include "Minmatching/PerfectMatching.h" 
#include <vector>
#pragma once

class MST {
public:
	float** adjacentMatrix;
	int* parent; //Array to store constructed MST
	int* key; //Key values used to pick minimum weight edge in cut
	bool* mstSet; //To represent set of vertices not yet included in MST
	int N; //the size of pointset
	//vector<vector<int>> TSP15;  

	MST(float** adjacentMatrix, int size);
	~MST();

	//deliverable a
	void makeTree();
	void printMST();

	//deliverable b
	float makeTSP2();
	stack<int> explore(int v, stack<int>); 

	//deliverable c
	float makeTSP1_5();
	
	float calMean(int option);
	float calStd(int option);

private:
	void minimumMatching();
	void combine(std::vector<int> oddDegree, PerfectMatching* pm, int node_num);
	int minKey(int key[], bool mstSet[]);

	void LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, vector<int> oddDegree, int N);

	void PrintMatching(int node_num, PerfectMatching* pm, vector<int> oddDegree);

	void printEulerUtil(int u, vector<vector<int>> adjList, bool visited[], stack<int> s); 
	bool isValidNextEdge(int u, int v, vector<vector<int>> adjList); 
	int DFSCount(int v, bool visited[], vector<vector<int>> adjLisit);

	void setVisited(int v); 
	bool isVisited(int v);  

};
