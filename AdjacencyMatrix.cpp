/***************************************************************
Title: AdjacencyMatrix.cpp
Author: Andres Quinones
Created on: December 13th, 2016
Description: AdjacencyMatrix for RNA Refragmentation Solution
Purpose: Graph Theory Independent Study
Usage: 
	to be used with the included main.cpp

Executable name is AllEuler unless desired otherwise.

***************************************************************/

/*
 * C++ Program to Implement Adjacency Matrix
   Inspired By Implimentation found on the following site
   http://www.sanfoundry.com/cpp-program-implement-adjacency-matrix/

   The online implimentation was inadequate for the purposes of this assignment
   It was heavily modified accordingly
	
	Euler circuits algorithm is within this class
 */
#include <iostream>
#include <vector>
#include <stack>
#include <string>
#include <cstdlib>
using namespace std;

// At first this struct was needed to store visited values for each edge
// But then the switch to using recursion made this no longer necessary
// Each adj[i][j] has a bool value for if it was visited and a string denoting the label of the edge 
struct Bool_String
{
	Bool_String(){};
	Bool_String& operator=(Bool_String &other){
		this->visited = other.visited;
		this->label = other.label;
		return *this;
	}
	bool visited;
	string label;
};

	// Edge is a simple class just for i,j,k triple 
	// (uniquely identifies edges and encapsulates all information needed)
class Edge
{
public:
	Edge(){};
	Edge(int x, int y, int z)
	{
		this->i = x;
		this->j = y;
		this->k = z;
	};
	void print()
	{
		cout << i << "-(" << k  << ")->" << j ;
	}

	int i;
	int j;
	int k;
};

/*
 * Adjacency Matrix Class
 */
class AdjacencyMatrix
{
    private:
    int n;
    vector<Bool_String> **adj;

    public:
	AdjacencyMatrix() {}
	// From a given n make 2d array of vectors of strings 
	//( this will represent our edges and labels for them)
    AdjacencyMatrix(int n)
    {
        this->n = n;
        adj = new vector<Bool_String>* [n];
        for (int i = 0; i < n; i++)
        {
            adj[i] = new vector<Bool_String> [n];
            for(int j = 0; j < n; j++)
            {
            	vector<Bool_String> x;
                adj[i][j] = x;
            }
        }
    }

    // Overloaded assignment. Copys all the data over this is instance
    AdjacencyMatrix operator=(AdjacencyMatrix &arg) // copy/move constructor is called to construct arg
	{
		this->n = arg.n;
		for (int i=0; i< this->n; i++)
			for (int j=0; j< this->n; j++)
			{
				this->adj[i][j] = arg.adj[i][j];
			}

		return *this;
	} // destructor of arg is called to release the resources formerly held by *this

    /*
     * Adding Edge to Graph
     */ 
    void add_edge(int origin, int destin, string s)
    {
        if( origin > n || destin > n || origin < 0 || destin < 0)
        {   
            cout<<"Invalid edge!\n";
        }
        else
        {
        	Bool_String x;
        	x.visited = false;
        	x.label = s;
            adj[origin][destin].push_back(x);
        }
    }
    /*
     * Print the graph
     */ 
    void display()
    {
        int i,j,k;
        for(i = 0;i < n;i++)
        {
            for(j = 0; j < n; j++)
            {
            	cout << i << "," << j << ": { ";
            	if (adj[i][j].size() == 0)
            		cout << " NE ";
            	else
            	{
            		for(k=0; k < adj[i][j].size(); k++)
            		{
                  	  cout << adj[i][j][k].label;
                  	  if (k != adj[i][j].size()-1)
                  	  	cout << " | ";
            		}
            	}
            	cout<< " }"<< endl;
            }
            cout <<endl;
        }
    }

    // Accessor
    string at(int i, int j, int k)
    {
    	return adj[i][j][k].label;
    }
    int countEdges()
    {
    	int count=0;
    	for (int i=0; i < n; i++)
    	{
    		for (int j=0; j < n; j++)
    			count += adj[i][j].size();
    	}
    	return count;
    }


	// Helper function to compute euler circuits
	void getAllEulerCircuits_Rec(int start,int end, int max, vector<Edge> pastEdges, vector<vector<Edge> > &output)
	{
		// If we traversed max edges in graph. Push pastEdges into output.
		if (start == end && pastEdges.size() == max)
		{
			output.push_back(pastEdges);
			return;
		}

		// Row in adj is all edges that leave from start
		for (int j=0; j<n; j++)
		{
			// If no edge from i to j. try next vertex
			if (adj[start][j].empty())
			{
				continue;
			}
			else
			{
				// Iterate through all the edges from i to j
				for (int k=0; k < adj[start][j].size(); k++)
				{
					// Edge is a simple class just for i,j,k triple 
					// (uniquely identifies edges and encapsulates all information needed)
					Edge e;
					e.i = start;
					e.j = j;
					e.k = k;
					bool inPath = 0;
					// Check if the edge k from i to j is a pastEdge 
					// If pastEdge empty. Dont loop
					for (int h = 0; h < pastEdges.size(); h++)
					{
						if (e.i == pastEdges[h].i &&
							e.j == pastEdges[h].j &&
							e.k == pastEdges[h].k)
						{
							inPath = 1;
							break;
						}
					}
					// If that edge was not yet traversed
					if (!inPath)
					{
						// Mark edge visited. ---- NO LONGER NEEDED ---
						adj[start][j][k].visited = 1;
						// push that edge on path
						pastEdges.push_back(e);

						// Recursive call with start at j with new pastEdges
						getAllEulerCircuits_Rec(j,end, max, pastEdges,output);

						// Remove the edge tried from pastEdges
						pastEdges.pop_back();
						//Mark edge unvisited---- NO LONGER NEEDED ---
						adj[start][j][k].visited = 0;
					}
				}
			}
		}
	}

// Compute all euler circuits
	vector<vector<Edge> > getAllEulerCircuits(int start)
	{
		// Count all edges
		int allEdgeCount = countEdges();
		// Container for output
		vector<vector<Edge> > out;
		// Past edges empty.
		vector<Edge> past_e;
		//Call to recursive function. Find circuits that start and end at start node.
		getAllEulerCircuits_Rec(start, start, allEdgeCount, past_e,out);
		// Return container filled with all circuits from start to start.
		return out;
	}
};


























// Debug adjacency list 
/* int main()
{
    int nodes, max_edges, origin, destin;
    cout<<"Enter number of nodes: ";
    cin>>nodes;
    AdjacencyMatrix am(nodes);
    max_edges = nodes * (nodes - 1);
    for (int i = 0; i < max_edges; i++)
    {
        cout<<"Enter edge (-1 -1 to exit): ";
        cin>>origin>>destin;
        if((origin == -1) && (destin == -1))
            break;
        am.add_edge(origin, destin, "Hello");
    }
    am.display();
    return 0;
} */
