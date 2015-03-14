#include "MST.h"
#include <stack>
//#include "Minmatching/PerfectMatching.h"
#include <vector>
#include <algorithm>
#include <cstring>

MST::MST(float** input, int size) {
	adjacentMatrix = input;
	key = new int[size];   
        mstSet = new bool[size];  
	parent = new int[size];

	N = size;
}

MST::~MST() {

}

//use Prim's algorithm or Kruskal algorithm. Copied from 'http://www.geeksforgeeks.org/greedy-algorithms-set-5-prims-minimum-spanning-tree-mst-2/'
void MST::makeTree() { 
     // Initialize all keys as INFINITE
     for (int i = 0; i < N; i++)
        key[i] = INT_MAX, mstSet[i] = false;
 
     // Always include first 1st vertex in MST.
     key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
     parent[0] = -1; // First node is always root of MST 
 
     // The MST will have V vertices
     for (int count = 0; count < N-1; count++)
     {
        // Pick thd minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet);
 
        // Add the picked vertex to the MST Set
        mstSet[u] = true;
 
        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < N; v++)
           // mstSet[v] is false for vertices not yet included in MST
           // Update the key only if adjacentMatrix[u][v] is smaller than key[v]
          if (adjacentMatrix[u][v] && mstSet[v] == false && adjacentMatrix[u][v] <  key[v])
             parent[v]  = u, key[v] = adjacentMatrix[u][v];
     }

	int cost = 0; 

	// MST Cost
	for (int i = 0; i < N; i++){
	  cost+=adjacentMatrix[i][parent[i]]; 
	} 

	std::cout << "MST Cost " << cost << endl; 
}

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int MST::minKey(int key[], bool mstSet[])
{
   // Initialize min value
   int min = INT_MAX, min_index;
 
   for (int v = 0; v < N; v++)
     if (mstSet[v] == false && key[v] < min)
         min = key[v], min_index = v;
 
   return min_index;
}

// A utility function to print the constructed MST stored in parent[]
void MST::printMST() {
	cout<<endl;
	cout<<"Minimum spanning tree from the adjacency matrix"<<endl;
	cout<<"Edge   Weight"<<endl;
	for (int i = 1; i < N; i++) {
		cout<<parent[i]<<" - "<<i<<"  "<<adjacentMatrix[i][parent[i]]<<endl;
	}
}

//calculate mean of all edges in the MST
float MST::calMean(int option) {
	float mean = 0.0;

	if(option == MST_1) {
		//calculate
	}else if(option == TSP2) {

	} else if(option == TSP1_5) {

	}

	return mean;
}

//calculate standard deviation of all edges in the MST
float MST::calStd(int option) {
	float std = 0.0;

	if(option == MST_1) {
		//calculate
	}else if(option == TSP2) {

	} else if(option == TSP1_5) {

	}

	return std;
}

void MST::makeTSP2() {
	//make a Eulerian tour by DFS
	std::stack<int> s; 

	for (int i = 0; i < N; i++){
	  mstSet[i] = false; 
	}

	for (int v = 0; v < N; v++){
	  if (mstSet[v] == false)
            s = explore(v, s); 
	}
	
	parent[0] = s.top(); 
	
	while (s.empty() == false)
	{
	  s.pop(); 
	}

 	printMST(); 

	int cost = 0; 

	//calculate heuristic TSP cost
	for (int l = 0; l < N; l++){
          cost+=adjacentMatrix[l][parent[l]]; 
        }

	std::cout << "TSP2 Cost " << cost << endl;  
}

stack<int> MST::explore(int v, stack<int> s){
  mstSet[v] = true;
  bool end = true;  
 
  for (int k = 1; k < N; k++){
    // Explore adjacent vertices 
    if (parent[k] == v){
      end = false; 
      // If the adjacent vertices has not been visited 
      if (mstSet[k] == false){
        if (s.empty() == false){
          parent[k] = s.top(); 
          s.pop(); 
        }
        s = explore(k,s);  
      }
    }
  }

  // If dead end vertex reached
  if (end == true)
  {
    s.push(v); 
  }

  return s; 
}

void MST::makeTSP1_5() {
	
	//construct minimum-weight-matching for the given MST
	minimumMatching();

	//make all edges has even degree by combining mimimum-weight matching and MST
	//combine();

	//calculate heuristic TSP cost
	
	//Fleury's Algorithm 
}

// Print Euler tour tour starting from vertex u 
void MST::printEulerUtil(int u, vector<vector<int>> adjList, bool visited[], stack<int> s)
{
  // Recurse for all the vertices adjacent to this vertex 
  for (int i = 0; i < adjList[u].size(); i++){

    int v = adjList[u][i]; // Adjacent vertex 

    // Check if u-v edge is valid next edge (i.e. doesn't burn a bridge)
    if (v != -1 && isValidNextEdge(u,v, adjList) && 
	visited[v] == false /*&& isVisited(v) == false*/){

      if (isVisited(v) == true){
	s.push(u);

	// remove edges
	adjList[u][i] = -1; 

        for (int j = 0; j < adjList[v].size(); j++){
          if (adjList[v][j] == u){
            adjList[v][j] = -1; 
            break; 
          }
        }
	
	setVisited(v); 
	printEulerUtil(v, adjList, visited, s); 
	visited[v] = true; 
	continue; 
      }

      if (s.empty() == false){
	cout << v << " --- " << s.top() << endl; 
	s.pop(); 

        adjList[u][i] = -1; // Remove u-v edge 

        // Remove v-u edge 
        for (int j = 0; j < adjList[v].size(); j++){
          if (adjList[v][j] == u){
            adjList[v][j] = -1; 
          break; 
          }
        }

	setVisited(v); 
	printEulerUtil(v, adjList, visited, s); 
	visited[v] = true; 
	continue; 
      }

      cout << u << " - " << v << endl;
      
      adjList[u][i] = -1; // Remove u-v edge 

      // Remove v-u edge 
      for (int j = 0; j < adjList[v].size(); j++){
        if (adjList[v][j] == u){
          adjList[v][j] = -1; 
          break; 
        }
      }

      setVisited(v); 
      
      printEulerUtil(v, adjList, visited, s); 
      visited[v] = true; 
    }
  }
}

void MST::setVisited(int v)
{
  mstSet[v] = true; 
}

bool MST::isVisited(int v) 
{
  return mstSet[v]; 
}

bool MST::isValidNextEdge(int u, int v, vector<vector<int>> adjList)
{
  // The edge u-v is valid in one of the following two cases: 

  // 1) If v is the only adjacent vertex of u 
  int count = 0; // To store count of adjcent vertices 

  for (int i = 0; i < adjList[u].size(); i++){
    if (adjList[u][i] != -1)
      count++; 
  }

  if (count == 1)
    return true; 

  // If there are multiple adjacents, then u-v is not a bridge
  // Do following steps to check if u-v is a bridge 

  // 2.a) count vertices reachable from u 
  bool visited[N]; 
  memset(visited, false, N); // This is being reset

  int count1 = DFSCount(u, visited, adjList);

  // 2.b) Remove edge (u,v) and after removing edge, count 
  // vertices reachable from u  
  int i;  
  int j; 

  for (i = 0; i < adjList[u].size(); i++){
    if (adjList[u][i] == v){
      adjList[u][i] = -1; 
      break; 
    }
  }

  for (j = 0; j < adjList[v].size(); j++){
    if (adjList[v][j] == u){
      adjList[v][j] = -1; 
      break; 
    }
  }

  memset(visited, false, N); 
  int count2 = DFSCount(u, visited, adjList); 

  // 2.c) Add the edge back to the graph 
  adjList[u][i] = v; 
  adjList[v][j] = u; 
  
  // 2.d) If count1 is greater, then edge (u,v) is a bridge 
  return (count1 > count2)? false: true; 
}

// A DFS based function to count reachable vertices from v
int MST::DFSCount(int v, bool visited[], vector<vector<int>> adjList)
{
  visited[v] = true; 
  int count = 1; 

  // Recurse for all vertices adjacent to this vertex 
  for (int i = 0; i < adjList[v].size(); ++i){
    if (adjList[v][i] != -1 && visited[adjList[v][i]] == false)
      count+=DFSCount(adjList[v][i], visited, adjList); 
  }
  return count; 
}

void MST::minimumMatching() { //if you choose O(n^2)
	int degree; 
	int numOdd = 0; 
	std::vector<int> oddDegree; 	 

	//find minimum-weight matching for the MST. 
	for (int i = 0; i < N; i++){ //Go through each vertex
	  degree = 1; 
	  for (int j = 1; j < N; j++){ //Compute degree
	    if (parent[j] == i){
	      degree++; 
	    }
          }
	  // Account for extra count on vertex zero 
	  if (i == 0)
	    degree--; 
	  // If degree of vertex is odd, add to set of odd degree vertices
	  if (degree & 1){
	    oddDegree.push_back(i); 
	    numOdd++; 
          }	
	}	

	//you should carefully choose a matching algorithm to optimize the TSP cost.
	struct PerfectMatching::Options options; 
	int i, e, node_num = numOdd, edge_num = numOdd*(numOdd-1)/2; 
std::cout << "node_num: " << node_num << endl; 
	int* edges; 
	int* weights; 
	PerfectMatching *pm = new PerfectMatching(node_num, edge_num); 

	LoadInput(node_num, edge_num, edges, weights, oddDegree, numOdd); 

	for(e = 0; e < edge_num; e++){
	  pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]); 
	}

	pm->options = options; 
	pm->Solve();

	double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm); 

	PrintMatching(node_num, pm, oddDegree); 

	printf("Total cost of the perfect min-weight matching = %.1f\n", cost);  

	combine(oddDegree, pm, node_num);  
}

void MST::LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, vector<int> oddDegree, int N){
	int e = 0; 
	node_num = N; 
	edge_num = N*(N-1)/2; 

	edges = new int[2*edge_num]; 
	weights = new int[edge_num]; 

	for (int i = 0; i < N; ++i){
	  for (int j = i+1; j < N; ++j){
	    edges[2*e] = i; //index of vertex of odd degree 
	    edges[2*e+1] = j; 
	    weights[e] = adjacentMatrix[oddDegree[i]][oddDegree[j]]; 
	    e++; 
	  }
	}
	
	if (e != edge_num){
	  cout << "the number of edge is wrong" << endl;
	  exit(1); 
	}
}

void MST::PrintMatching(int node_num, PerfectMatching* pm, vector<int> oddDegree){
	int i,j; 
	
	for (i = 0; i < node_num; i++){
	  j = pm->GetMatch(i); 
	  if (i < j) printf("%d %d\n", oddDegree[i], oddDegree[j]); 
	}
}

void MST::combine(std::vector<int> oddDegree, PerfectMatching* pm, int node_num) {
	// Data struture to hold combined MST and minimum matching edges
	std::vector<std::vector<int>> adjList(N); 
	int j; // Contains perfect matching vertex 

	//Fill in new data structure with MST 
	for (int i = 1; i < N; i++){
	  adjList[i].push_back(parent[i]); 
	  adjList[parent[i]].push_back(i); 
	}

	//combine minimum-weight matching with the MST to get a multigraph which has vertices with even degree
	for (int i = 0; i < node_num; i++){
	  j = pm->GetMatch(i); 
	  if (i < j){
	    adjList[oddDegree[j]].push_back(oddDegree[i]); 	
	    adjList[oddDegree[i]].push_back(oddDegree[j]); 
	  }
	}

 
	std::cout << "Printing Eulerian circuit" << endl; 

/*
std::cout << "Adjacent List size: " << adjList.size() << endl; 

	// Print out to see if it is accurate
	for (int i = 0; i < N; i++){
	  for (int j = 0; j < adjList[i].size(); j++){
	    //cout << "Number of vertices" << adjList[i].size(); 
	    cout << i << " - " << adjList[i][j] << endl; 
	  }
	} 
*/
	bool visited[N]; 
 	memset(visited, false, N);

	memset(mstSet, false, N);  

	stack<int> s; 

	printEulerUtil(0, adjList, visited, s); 

}
