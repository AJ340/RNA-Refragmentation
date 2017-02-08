/***************************************************************
Title: main.cpp
Author: Andres Quinones
Created on: December 13th, 2016
Description: RNA-Refragmentation Solution / Demo
Purpose: Graph Theory Independent Study
Usage: 
	make
								or
	g++ -std=c++11 main.cpp -o RNA-Refragmentation-AQ

Executable name is RNA-Refragmentation-AQ unless desired otherwise.

***************************************************************/




#include "RNA-StringParse.cpp"
#include <unordered_map>
#include "AdjacencyMatrix.cpp"

// FUNCTION DECLARATIONS START. DEFINITIONS BELOW MAIN
//
// Builds unordered map from non_single.
// So that we can lookup a index in adjacency matrix by the string.
vector<string> buildVertices(	unordered_map <string,int> &vertexes, vector <string> non_single);
//
// Build ajacency matrix from non_single
// The size of the unordered map is the size of the adj matrix.
void build_AM_From_Non_Single(AdjacencyMatrix &am, const unordered_map <string,int> vertexes, vector <string> non_single);
//
// Adjacency Matrix returns i,j,k values for edges in the graph. 
// They must be decoded to the strings represented by the vertexes and edges
vector<string> decodeEulerCircuits(vector<vector<Edge> > eulerCircuits, AdjacencyMatrix a, vector<string> ind_to_string);
//
// Helper to remove x's and /
// x represents No Edge in graph
bool IsSlashOrX(char c);
//
// Removes slashes and x's from a string
void RemoveSlashX(string &x);
//
// Get rid of all edges that dont end in the computed end value
// Also removes the / and x's from all of the circuits
vector<string> cleanEulerCircuits (vector<string> &input, string end);
//
// This function is where the examples are set. 
vector<string> init_Examples( );
//
// This is a DEBUG function that tests just the example from the notes.
int test_Notes_Example(int argc, char * argv[]);
//
// FUNCTION DECLARATIONS END.






// Main function. Run on all examples set by init_Examples()
int main(int argc, char * argv[])
{
	// Get examples
	vector<string> input_examples = init_Examples();

	for( int i=0; i < input_examples.size(); i+=2)
	{
		// All of our needed variables
		string x = input_examples[i];
		string y = input_examples[i+1];
		string start;
		string end;
		vector<string> ind_to_string;
		unordered_map<string,int> am_vertices;
		vector<vector<Edge> > eulerCircuits;
		vector<string> decodedCircuits;
		vector<string> cleanDecodedCircuits;

		cout << "--------------------------- Start of Example " << i/2+1 << " ------------------------" << endl << endl;
		cout << "Target RNA After G_Enzyme:  " << x << endl;
		cout << "Target RNA After UC_Enzyme:  " << y << endl;

		cout << "Ammount of non_single fragments: ";
		vector <string> non_single = processInputEnzymes(x, y, start, end);
		cout << non_single.size() << endl;

		cout << "\n Ammount of vertexes made from non_single fragments: ";
		ind_to_string = buildVertices(am_vertices,non_single);
		cout << am_vertices.size() << endl;
		AdjacencyMatrix am(am_vertices.size());

		build_AM_From_Non_Single(am, am_vertices, non_single);

		cout << "RNA reconstruction must begin with [" << start << "] and end with [" << end << "] \n"; 

		cout << "_____________________________________________"<< endl;
		cout << "||Index |  Vertex Label.                     " << endl  ;
		for(int i = 0; i < ind_to_string.size(); i++)
		{
			cout << "|| "<< i << "    | " << ind_to_string[i] << endl;
		}
		cout << "_____________________________________________"<< endl << endl;
		cout << "Adjacency Matrix Displayed: \n";
		am.display();
		cout << endl << endl;


		int start_index = am_vertices[start];
		eulerCircuits = am.getAllEulerCircuits(start_index);
		decodedCircuits = decodeEulerCircuits(eulerCircuits, am, ind_to_string);
		cleanDecodedCircuits = cleanEulerCircuits(decodedCircuits, end);

		cout << "All Euler Circuits Found: \n";
		for (int i = 0; i < cleanDecodedCircuits.size(); i++)
			cout << cleanDecodedCircuits[i] << endl;
		cout << "Total Found " << cleanDecodedCircuits.size() << endl<< endl;
		cout << "--------------------------- End of Example " << i/2+1 << " ------------------------" << endl << endl;


	}

}

// Builds unordered map from non_single.
// So that we can lookup a index in adjacency matrix by the string.
vector<string> buildVertices(	unordered_map <string,int> &vertexes, vector <string> non_single)
{
	// Parse non_single and make vertexes for the first and last of each fragment
	vector<string> s;
	int index = 0;
	for (int i=0; i < non_single.size(); i++)
	{
		string previous = "", token;
		istringstream ss(non_single[i]);
		getline(ss, token, '/');
		// Make vertex from first
		if (vertexes.find(token) == vertexes.end())
		{
			vertexes.insert({token, index});
			s.push_back(token);
			index++;
		}
		//cout << "First: " << first;
		while(getline(ss, token, '/'))
		{
			previous = token;
		}

		if (vertexes.find(token) == vertexes.end())
		{
			vertexes.insert({previous,index});
			s.push_back(token);
			index++;
		}
	}
	return s;
}


// Build ajacency matrix from non_single
// The size of the unordered map is the size of the adj matrix.
void build_AM_From_Non_Single(AdjacencyMatrix &am, const unordered_map <string,int> vertexes, vector <string> non_single)
{
	int index = 1, i=0, j=0;

	unordered_map<string,int>::const_iterator c_i;
	// For each fragment..
	for (int w=0; w < non_single.size(); w++)
	{
		string edges = "x";
		string previous = "", token;
		istringstream ss(non_single[w]);
		getline(ss, token, '/');
		// Find from vertex for that edge
		i = vertexes.at(token);

		//cout << "First: " << first;
		while(getline(ss, token, '/'))
		{
			if (previous != "")
				edges += previous+"/";
				//cout << " Inner: " << previous;
			previous = token;
		}
		// Find to vertex for that edge
		j = vertexes.at(previous);
		
		if (edges.size() == 1)
			// Edges is "x". This will denote an edge without a label
			am.add_edge(i,j,edges);
		else
			// Otherwise label the edge the accumulated string except cut off the extra / at the end
			am.add_edge(i,j,edges.substr(1,edges.size()-2));

		//Make vertex from last
		//cout << " End: " << previous << endl;
	}
	return;
}

// Adjacency Matrix returns i,j,k values for edges in the graph. 
// They must be decoded to the strings represented by the vertexes and edges
vector<string> decodeEulerCircuits(vector<vector<Edge> > eulerCircuits, AdjacencyMatrix a, vector<string> ind_to_string)
{
	vector<string> decoded;

	for (int i = 0; i < eulerCircuits.size(); i++)
	{
		string s;
		int beggining = 0;
		int ec_size = eulerCircuits[i].size();
		for (int j = 0; j < eulerCircuits[i].size(); j++)
		{
			if (beggining == 0)
				s += ind_to_string [eulerCircuits[i][j].i];
			s += a.at(eulerCircuits[i][j].i,eulerCircuits[i][j].j,eulerCircuits[i][j].k);
			if (beggining != ec_size-1)
			s += ind_to_string[eulerCircuits[i][j].j];
			beggining ++;
		}
		decoded.push_back(s);
	}
	return decoded;
}

// Helper to remove x's and /
// x represents No Edge in graph
bool IsSlashOrX(char c)
{
    switch(c)
    {
    case '/':
    return true;
    case 'x':
        return true;
    default:
        return false;
    }
}

// Removes slashes and x's from a string
void RemoveSlashX(string &x)
{
	string s;
	for (const auto c: x)
		if (!IsSlashOrX(c))
			s.push_back(c);
	x = s;
}

// Get rid of all edges that dont end in the computed end value
// Also removes the / and x's from all of the circuits
vector<string> cleanEulerCircuits (vector<string> &input, string end)
{
	vector<string> clean;
	RemoveSlashX(end);
	for (int i = 0; i < input.size(); i++)
	{
		RemoveSlashX(input[i]);
		string in_end = input[i].substr(input[i].size()-end.size(), input[i].size());
		if (in_end == end)
			clean.push_back(input[i]);
	}
	return clean;
}

// This function is where the examples are set. 
vector<string> init_Examples( )
{
	// For each example we push G first and then UC after
	vector<string> examples;

	// NOTES EXAMPLE
	examples.push_back("AUCG,AUG,G,CU,ACUAUACG");
	examples.push_back("GGAC,U,AU,GAU,C,U,AC,GC,AU");
	//

	// EXAMPLE 1
	// Single Strand PR0021. According to http://ndbserver.rutgers.edu/
	// STRUCTURE OF THE RIBONUCLEOPROTEIN CORE OF THE E. COLI SIGNAL RECOGNITION PARTICLE
	string self_example1 = "GGCUCUGUUUACCAGGUCAGGUCCGAAAGGAAGCAGCCAAGGCAGAGCCCC";
	G_Enzyme(self_example1,',');
	examples.push_back(self_example1);
	//
	//
	self_example1 = "GGCUCUGUUUACCAGGUCAGGUCCGAAAGGAAGCAGCCAAGGCAGAGCCCC";
	UC_Enzyme(self_example1,',');
	examples.push_back(self_example1);
	//

	// EXAMPLE 2
	// According to http://ndbserver.rutgers.edu/. 5TPY. 
	// EXONUCLEASE RESISTANT RNA FROM ZIKA VIRUS
	// NOTE : THIS EXAMPLE NOT USED BECAUSE IT WAS TOO COMPLEX. GRAPH CREATED IS GIVEN THOUGH :)
	//
	string self_example2 = "GGGUCAGGCCGGCGAAAGUCGCCACAGUUUGGGGAAAGCUGUGCAGCCUGUAACCCCCCCACGAAAGUGGG";
	G_Enzyme(self_example2,',');
	examples.push_back(self_example2);
	//
	self_example2 = "GGGUCAGGCCGGCGAAAGUCGCCACAGUUUGGGGAAAGCUGUGCAGCCUGUAACCCCCCCACGAAAGUGGG";
	UC_Enzyme(self_example2,',');
	examples.push_back(self_example2);
	//
	//

	// EXAMPLE 3
	// This example is made up.
	string self_example3 = "AAAUCAUAGGCAGAAUGCAUUCGAUGCAAU";
	G_Enzyme(self_example3,',');
	examples.push_back(self_example3);
	//
	self_example3 = "AAAUCAUAGGCAGAAUGCAUUCGAUGCAAU";
	UC_Enzyme(self_example3,',');
	examples.push_back(self_example3);
	//


	return examples;
}







// This is a DEBUG function that tests just the example from the notes.
int test_Notes_Example(int argc, char * argv[])
{
	string x = "AUCG,AUG,G,CU,ACUAUACG";
	string y = "GGAC,U,AU,GAU,C,U,AC,GC,AU";
	string start;
	string end;
	vector<string> ind_to_string;
	unordered_map<string,int> am_vertices;
	vector<vector<Edge> > eulerCircuits;
	vector<string> decodedCircuits;

	vector <string> non_single = processInputEnzymes(x, y, start, end);
	cout << endl << non_single.size() << endl;

	ind_to_string = buildVertices(am_vertices,non_single);
	cout << endl << am_vertices.size() << endl;
	AdjacencyMatrix am(am_vertices.size());
	AdjacencyMatrix am_copy(am_vertices.size());


	build_AM_From_Non_Single(am, am_vertices, non_single);
		am_copy = am;
	cout << "AM COPY DISPLAY";
	am_copy.display();
	cout << "NDIOBUFWF " << start << end << endl;
	int start_index = am_vertices[start];
	int end_index = am_vertices[end];
	cout << "START AND END INDEX SUCCESS " << start << end << endl;
	am.display();
	cout << " Vertex :" << am_vertices.at("AU") << " " << am_vertices.at("C") << "\n";

	eulerCircuits = am.getAllEulerCircuits(0);
	decodedCircuits = decodeEulerCircuits(eulerCircuits, am, ind_to_string);
	decodedCircuits = cleanEulerCircuits(decodedCircuits, end);
	cout << "All Euler Circuits Found: ";
	for (int i = 0; i < decodedCircuits.size(); i++)
		cout << decodedCircuits[i] << endl;
	return 0;
}