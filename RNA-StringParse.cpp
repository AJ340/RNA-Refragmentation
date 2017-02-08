/***************************************************************
Title: RNA-StringParse.cpp
Author: Andres Quinones
Created on: December 13th, 2016
Description: Parses RNA Fragment Strings
Purpose: Graph Theory Independent Study
Usage: 
	to be used with the included main.cpp

Executable name is AllEuler unless desired otherwise.

***************************************************************/

#include "RNA-Utility.cpp"

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

// type = what enzyme was just run on input
// If any end in a character the enzyme breaks, its abnormal (it must be end)
string findAbnormal(string type, string input) 
{
	istringstream ss1(input);
	string token, empty;
	while(getline(ss1, token, ','))
	{
		if (token.empty())
			return empty;
		else if (type == "G" && token.back() != 'G')
			return token;
		else if (type == "UC" && (token.back() != 'U' && token.back() != 'C'))
			return token;
	}
	return empty;
}

// Given the inners and singles. Returns the ones that arent in both
// Expected to be 2 if used as intended in notes
vector<string> computeStartAndEnd(vector<string> inner, vector<string> single)
{
	vector<string> se;
	string start;
	int i, j, s_size = single.size(), i_size = inner.size();
	for (i=0; i < s_size; i++)
	{
		bool found = 0;
		for (j=0; j < i_size;j++)
		{
			if (single[i] == inner[j])
			{
				found = 1;
				inner.erase(inner.begin()+j);
				single.erase(single.begin()+i);
				i_size--;s_size--;i--;j--;
				break;
			}
		}
		if (!found)
		{
			se.push_back(single[i]);
		}
	}
	return se;
}

// Single and Inners are returned in params
void getSingleAndInner(string input, vector<string> &inner_bases, vector<string> &single_bases, vector<string> &some_nonsingle) 
{
	istringstream ss1(input);
	string token1, token2, first;
	while(getline(ss1, token1, ','))
	{
		int fragment_count = 1;
		string previous = "";
		istringstream ss2(token1);
		getline(ss2, first, '/');
		// Make vertex from first

		//cout << "First: " << first;
		while(getline(ss2, token2, '/'))
		{
			fragment_count++;
			if (previous != "")
				inner_bases.push_back(previous);
				//cout << " Inner: " << previous;
			previous = token2;
		}
		if (fragment_count == 1)
			single_bases.push_back(first);
		else 
			some_nonsingle.push_back(token1);
		//Make vertex from last
		//cout << " End: " << previous << endl;
	}
}

// Expected possible starts to be 2
// Erases the one that matches our computed end
// The one left is the start string
string findStartFrom_SE (vector<string> possible_starts, string end)
{
	int index;
	vector<string> empty = possible_starts;
	for (int i=0; i < possible_starts.size(); i++)
	{
		index = end.find(possible_starts[i]);
		if (index != string::npos)
		{
			empty.erase(empty.begin()+i);
			break;
		}
	}
	return empty[0];
}

// Adds Start after End in nonsingles
void enlargeEnd ( vector<string> &non_single, string start, string end)
{
	for (int i = 0; i< non_single.size(); i++)
	{
		if (non_single[i] == end)
			non_single[i] = end + "/" + start;
	}
}

// X = fragments after G enzyme
// Y = fragments after UC enzyme
// Computes start_node and end node
// Returns the non single fragments
vector<string> processInputEnzymes(string x, string y, string &start_node, string &end_node)
{
	vector<string> inner, single, s_and_e, non_single;

	end_node = findAbnormal("G",x);
	if (end_node.empty())
		end_node = findAbnormal("UC",y);
	if (end_node.empty())
		cout << "Error \n";

    UC_Enzyme(x);
    G_Enzyme(y);

    // Compute Single and Inner on both
    getSingleAndInner(x,inner,single, non_single);
    getSingleAndInner(y,inner,single, non_single);
	
    s_and_e = computeStartAndEnd(inner, single);

    start_node = findStartFrom_SE(s_and_e, end_node);
	
	UC_Enzyme(end_node);
	G_Enzyme(end_node);

	enlargeEnd(non_single, start_node, end_node);
	
	return non_single;
}




// This is a debug function to ensure everything included works
void test_String_Parse(string x, string y)
{
	vector<string> inner, single, s_and_e, non_single;
	string end_node, start_node;

	end_node = findAbnormal("G",x);
	if (end_node.empty())
		end_node = findAbnormal("UC",y);
	if (end_node.empty())
		cout << "Error \n";

	cout << "\n End_Node: "<< end_node << endl;
    UC_Enzyme(x);
    G_Enzyme(y);

    cout << "\n G-Enzyme input after UC_Enzyme: \n" << x << endl;
    cout << "\n UC-Enzyme input after G_Enzyme: \n" << y << endl;


    // Compute Single and Inner on both
    getSingleAndInner(x,inner,single, non_single);
    getSingleAndInner(y,inner,single, non_single);
	cout << "\n Inner Bases: ";
    for (int i=0; i < inner.size(); i++)
    	cout << inner[i] << " ";
    cout << endl;

    cout << "\n Single Bases: ";
    for (int i=0; i < single.size(); i++)
    	cout << single[i] << " ";
    cout << endl;

    cout << "\n Non-Single Bases: ";
    for (int i=0; i < non_single.size(); i++)
    	cout << non_single[i] << " ";
    cout << endl;



    s_and_e = computeStartAndEnd(inner, single);
    cout << "\n Start and End: ";
    for (int i=0; i < s_and_e.size(); i++)
    	cout <<  s_and_e[i] << " ";
    cout << endl;


    start_node = findStartFrom_SE(s_and_e, end_node);
	cout << "Start Node: "<< start_node << endl;

	UC_Enzyme(end_node);
	G_Enzyme(end_node);

	cout << "End: " << end_node << endl;

	enlargeEnd(non_single, start_node, end_node);
	cout << "\n Enlarged end: ";
    for (int i=0; i < non_single.size(); i++)
    	cout <<  non_single[i] << " ";
    cout << endl;
	return;
}