#include <fstream>
#include <iostream>
using namespace std;

void simc_replay_randomseed(string infile, string outfile, string word, string wordreplace)
{
	ifstream in(infile);
	ofstream out(outfile);
	string wordToReplace(word);
	string wordToReplaceWith(wordreplace);

	if (!in)
	{
		cerr << "Could not open " << infile << "\n";
		return 1;
	}

	if (!out)
	{
		cerr << "Could not open " << outfile << "\n";
		return 1;
	}

	string line;
	size_t len = wordToReplace.length();
	while (getline(in, line))
	{
		while (true)
		{
			size_t pos = line.find(wordToReplace);
			if (pos != string::npos)
				line.replace(pos, len, wordToReplaceWith);
			else 
				break;
		}

		out << line << '\n';
	}
}
