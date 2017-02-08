/***************************************************************
Title: RNA-Utility.cpp
Author: Andres Quinones
Created on: December 13th, 2016
Description: G and UC Enzyme Simulation
Purpose: Graph Theory Independent Study
Usage: 
	to be used with the included main.cpp

Executable name is AllEuler unless desired otherwise.

***************************************************************/

#include <iostream>
#include <string>
using namespace std;

// Runs UC_Enzyme seperates splits by '/'
void UC_Enzyme(string &G_Frag)
{
	for (int i=0; i<G_Frag.length(); i++)
	{
		if(G_Frag[i] == 'U' || G_Frag[i] == 'C')
		{
			if (i < G_Frag.length()-1 && G_Frag[i+1] != ',')
			{
				G_Frag.insert(i+1, "/");
				i++;
			}
		}
	}
	return;
}

// Runs UC_Enzyme seperates splits by seperator
void UC_Enzyme(string &G_Frag, char seperator)
{
	for (int i=0; i<G_Frag.length(); i++)
	{
		if(G_Frag[i] == 'U' || G_Frag[i] == 'C')
		{
			if (i < G_Frag.length()-1 && G_Frag[i+1] != ',')
			{
				string s;
				s.push_back(seperator);
				G_Frag.insert(i+1, s);
				i++;
			}
		}
	}
	return;
}

// Runs G_Enzyme seperates splits by '/'
void G_Enzyme(string &UC_Frag)
{
	for (int i=0; i<UC_Frag.length(); i++)
	{
		if(UC_Frag[i] == 'G')
		{
			if (i < UC_Frag.length()-1 && UC_Frag[i+1] != ',')
			{
				UC_Frag.insert(i+1, "/");
				i++;
			}
		}
	}
	return;
}	

// Runs UC_Enzyme seperates splits by seperator
void G_Enzyme(string &UC_Frag, char seperator)
{
	for (int i=0; i<UC_Frag.length(); i++)
	{
		if(UC_Frag[i] == 'G')
		{
			if (i < UC_Frag.length()-1 && UC_Frag[i+1] != ',')
			{
				string s;
				s.push_back(seperator);
				UC_Frag.insert(i+1, s);
				i++;
			}
		}
	}
	return;
}	