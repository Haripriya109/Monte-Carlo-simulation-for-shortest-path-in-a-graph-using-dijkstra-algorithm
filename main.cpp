#include<bits/stdc++.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <vector>

using namespace std;

const int num_of_vertices = 7;
const int inf = 999;

//Graph Class
//random_graph method generates the graph.
class Graph {
public:
    void random_graph(float edge_density, int distance_range[]);
    void add_edge (int node_1, int node_2);
    bool adjacent (int node_1, int node_2);
    vector<int> neighbours (int node_1);
    void delete_edge (int node_1, int node_2);
    int get_node_value (int node_1);
    void set_node_value (int node_val, int node_1);
    int get_edge_value (int node_1, int node_2);
    void set_edge_value (int val, int node_1, int node_2);
    float prob_fn();
    void print_graph();

    vector<int> adj[num_of_vertices];
    int node_wt[num_of_vertices];
    int edge_wt[num_of_vertices][num_of_vertices];
};

//To add an edge to the graph.
//Inputs are the two nodes the new edge should be connected to.

void Graph::add_edge (int node_1, int node_2) {
    bool edge_exists = adjacent(node_1, node_2); //Check if the edge already exists.
    if (edge_exists == false) {
        adj[node_1].push_back(node_2);
        adj[node_2].push_back(node_1);
    }
}
//to delete an edge
//inputs are the two nodes that the edge is connected to.
void Graph::delete_edge(int node_1, int node_2) {
    bool edge_exists = adjacent(node_1, node_2); //Check if the edge exists
    if (edge_exists == true) {
        adj[node_1].erase (adj[node_1].begin()+node_2);
        adj[node_2].erase (adj[node_2].begin()+node_1);
    }
}
//Check if there is an edge between the two nodes, node_1 and node_2
bool Graph::adjacent (int node_1, int node_2) {
    if (find(adj[node_1].begin(), adj[node_1].end(), node_2) == adj[node_1].end()) {return false;}
    return true;
}
//Returns all of the nodes that node_1 is connected to, directly.
vector<int> Graph::neighbours (int node_1) {
    return adj[node_1];
}
//Returns the weight of node_1
int Graph::get_node_value (int node_1) {
    int retVal = node_wt[node_1];
    return retVal;
}
//To set the weight of the node
void Graph::set_node_value (int node_val, int node_1) {
    node_wt[node_1] = node_val;
}
//Returns the weight of the edge
int Graph::get_edge_value (int node_1, int node_2) {
    bool edge_exists = adjacent(node_1, node_2); //Check if the edge exisits
    if (edge_exists == true) {
        int retVal = edge_wt[node_1][node_2];
        return retVal;
    }
    return 0;
}
//Set the weight of the edge
void Graph::set_edge_value (int val, int node_1, int node_2) {
    bool edge_exists = adjacent(node_1, node_2);
    if (edge_exists == true) {
        edge_wt[node_1][node_2] = val;
        edge_wt[node_2][node_1] = val;
    }
}
//Generates Random Graph
void Graph::random_graph (float edge_density, int distance_range[]) {

    // To generate a random number
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distr(distance_range[0], distance_range[1]); //The range is from 1 - 10

    float prob_val = 0;

    for(int i= 0; i<num_of_vertices; i++) {
        for(int j = 0; j<num_of_vertices; j++) {
            prob_val = prob_fn();
            if (i != j && prob_val < edge_density) {  //Adding edge according to the density
                add_edge(i, j);
                set_edge_value(distr(gen), i, j); //Weight of the edge is randomly generated
            }
        }
    }
}
//Probability function
float Graph::prob_fn () {
  float pdf = rand() % 400;
  if (pdf > 360)
     return 1;
  else if (pdf < 0)
     return 0;
  else
     return pdf/360;
}
//To print the generated graph
void Graph::print_graph() {
	for (int v = 0; v < num_of_vertices; ++v) {
		cout << "\n Adjacency list of vertex " << v << "\n head ";
		for (auto x : adj[v]) {
            cout << "-> " << x;
		}
		cout<<"\n";
		//To print edge weight, uncomment the following block:
		//------------------------------------------------------------------------------------------------
		/*
		for (auto x : adj[v]) {
            cout << "Edge weight: (" << v << " to " << x << ") " << edge_wt[v][x] << "\n";
		}
		*/
		//-------------------------------------------------------------------------------------------------
	}
}

class ShortestPath {
public:
    int find_shortest_path (int src, int destn, Graph my_graph);
    void printShortestPath (vector<int> closed_set);
	vector<int> getShortestPath(vector<int> closed_set, Graph my_graph);
	int getIndex(vector<int> compVal, int minVal);
	int getPathCost(int node_wt[], int dest);
};
//To calculate the cost of the path
int ShortestPath::getPathCost(int node_wt[], int dest) {
    int retVal = node_wt[dest];
    if (node_wt[dest] != inf) {
        cout<<"\nCost found for the path is " << node_wt[dest];
        return retVal;
    }
    cout << "\nNo Path Found!" << endl;
    return 0;
}

//To get the shortest path from the closed set
vector<int> ShortestPath::getShortestPath(vector<int> closed_set, Graph my_graph) {

    vector<int> neighbours;
    vector<int> retSet;
    auto set_length  = closed_set.size();
    //cout << "\nset_length: " << set_length;
    int i = set_length - 1;
    retSet.push_back(closed_set[i]);
    auto itr = closed_set.end()-2;
    while (i != 0) {
        if(my_graph.adjacent(closed_set[i], *itr)) {
            neighbours = my_graph.neighbours(closed_set[i]);
            auto neighbours_len = neighbours.size();
            for(int j = 0; j < neighbours_len; j++) {
                if (find(retSet.begin(), retSet.end(), neighbours[j]) == retSet.end()) {
                    if (my_graph.node_wt[neighbours[j]]+my_graph.edge_wt[neighbours[j]][closed_set[i]] == my_graph.node_wt[closed_set[i]]) {
                        retSet.push_back(neighbours[j]);
                        i = getIndex(closed_set, neighbours[j]);
                        itr--;
                    }
                }
            }

        } else {
            itr--;
            continue;
        }
    }

    return retSet;
}

int ShortestPath::getIndex(vector<int> compVal, int minVal) {
    for (int i = 0; i < num_of_vertices; i++) {
        if (compVal[i] == minVal) {
            return i;
        }
    }
    return -1;
}
//This function returns the cost of the shortest path after finding the path.
//To debug, uncomment the print statements.
int ShortestPath::find_shortest_path(int src, int destn, Graph my_graph) {
    for (int i = 0; i < num_of_vertices; i++) {
        my_graph.node_wt[i] = inf;
    }
    my_graph.node_wt[src] = 0;

    vector<int> open_set, closed_set;
    closed_set.push_back(src);

	int count = 0;
	int node = src;
	int dist = 0;
	vector<int> compVal(num_of_vertices, inf);
	while(count  < num_of_vertices) {
		vector<int>::iterator it;
		//cout<<"\nAt node: " << node;
		for (it = my_graph.adj[node].begin(); it != my_graph.adj[node].end(); it++) {
            //cout << "\nChecking for node " << *it << "--->" << node;
			if (find(closed_set.begin(), closed_set.end(), *it) == closed_set.end()) {
				if (find(open_set.begin(), open_set.end(), *it) == open_set.end()) {
					open_set.push_back(*it);
					/*
					cout << "\nPrinting Open Set: ";
					for (auto x : open_set) {
                        cout<< "\n" << x;
					}
					*/
				}
				dist = my_graph.node_wt[node] + my_graph.edge_wt[*it][node];
				//cout << "\nDistance calculated: " << dist;
				if (dist < my_graph.node_wt[*it]) {
                    my_graph.node_wt[*it] = dist;
                    compVal[*it] = dist;
                    /*
                    cout << "\nPrinting node weight: ";
                    for (auto x : my_graph.node_wt) {
                        cout << "\n" << x;
                    }

                    cout << "\nPrinting CompVal mat: ";
                    for (auto x: compVal) {
                        cout << "\n" << x;
                    }
                    */
                }
			}
		}
		int minVal = *min_element(compVal.begin(), compVal.end());
		//cout << "\nMin val Calculated is : " << minVal;

		int node_val = getIndex(compVal, minVal);
        //cout << "\nMin val is at index: " << node_val;

		if (find(closed_set.begin(), closed_set.end(), node_val) == closed_set.end()) {
			closed_set.push_back(node_val);
			compVal[node_val] = inf;
			/*
			cout << "\nPrinting closed set: ";
			for (auto x : closed_set) {
                cout << "\n" << x;
			}
			cout << "\nPrinting comp Val mat: ";
			for (auto x: compVal) {
                cout << "\n" << x;
			}
			cout << "\n";
			cout << node << " | " << node_val << " | " << destn;
			*/
			if (find(closed_set.begin(), closed_set.end(), destn) != closed_set.end()) {
                //cout << "\nDestination " << destn << " Reached!";
				break;
			}
		}
		node = node_val;
		count ++;
	}
    int x = getPathCost(my_graph.node_wt, destn);
    if (x != 0) {
        vector<int> shortestPath = getShortestPath(closed_set, my_graph);
        //cout << "\nThe Shortest Path is: \n";
        cout << "\nThe shortest path found for source " << src << " and destination " << destn << " is: ";
        printShortestPath(shortestPath); //prints the shortest path.
    }

    return x;
}

void ShortestPath::printShortestPath (vector<int> shortestPath) {
    int i = shortestPath.size()-1;
    while ( i >= 0) {
        cout<< shortestPath[i] << "-->";
        i--;
    }
    cout << "end" << endl;

}

class MonteCarlo {
public:
    int runShortestPathAlgo (Graph my_graph, ShortestPath shortestPathObj, int src, int dest);
    float getAvgCost(int pathCostArr[]);
    void run();
};

//Runs the functions that calculates shortest path and its cost
int MonteCarlo::runShortestPathAlgo (Graph my_graph, ShortestPath shortestPathObj, int src, int dest) {
    int cost = shortestPathObj.find_shortest_path(src, dest, my_graph);

    return cost;
}
//This function calculates the average cost
float MonteCarlo::getAvgCost(int pathCostArr[]) {
    float sum = 0;
    for(int i = 0; i < num_of_vertices -1 ; i++) {
        sum = sum + pathCostArr[i];
    }
    float retVal = sum/(num_of_vertices-1);
    return retVal;
}

void MonteCarlo::run() {
    int src = 0; //Source is set to be 0
    float edge_densities[2] = {0.20, 0.40};
    int distance_range[2] = {1, 10};
    Graph my_graph;
    ShortestPath shortestPathObj;
    int pathCostArr[num_of_vertices-1];
    for (int d = 0; d < 2; d++) {
        int turn = d+1;
        cout << "\n-------------------------Running Simulation " << turn << " for density: " << edge_densities[d] << " ---------------------------" <<endl;
        my_graph.random_graph(edge_densities[d], distance_range);
        cout << "\nThe random graph generated for Dijkstra's shortest path algorithm is: "<< endl;
        my_graph.print_graph();
        for(int i = src + 1; i < num_of_vertices; i++) { //Considering source to be 0, calculate path from source to every node, one at a time.
            pathCostArr[i-1] = runShortestPathAlgo(my_graph, shortestPathObj, src, i);
        }
        float avgCost = getAvgCost(pathCostArr);
        cout << "\nAverage Cost is: " << avgCost << endl;
    }
}


int main () {
    MonteCarlo m;
    m.run();
}
