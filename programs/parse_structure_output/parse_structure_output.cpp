//Author: Sarah P. Flanagan
//Date: 10 March 2016
//Purpose: To take the output from structure (in the /Results/ directory) 
//and extract the info about the clusters


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <stdlib.h>

using namespace std;

bool FileTest(ifstream& file, string filename)
{
	cout << filename;
	if (file.is_open())
		cout << " open\n";
	else
	{
		while (!file.is_open())
		{
			cout << " not open. Please re-enter filename: ";
			getline(cin, filename, '\n');
			file.open(filename);
		}
	}
	return true;
}

istream& universal_getline(istream& is, string& t)
{
	//this code is adapted from a post on stackoverflow:
	// http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	//written by user763305
	t.clear();
	istream::sentry se(is, true);
	streambuf* sb = is.rdbuf();//sets pointer to stream buffer object

	for (;;)
	{
		int c = sb->sbumpc();//get current character and advance to the next position
		switch (c)//tests for equality against a list of variables (like multiple if statements)
		{
		case '\n'://if the next character is '\n', return the line
			return is;
		case '\r'://if the character is '\n', see if the next one is '\n' or return the line
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:// Also handle the case when the last line has no line ending
			if (t.empty())//if there's nothing there, set it to be the end of file
				is.setstate(ios::eofbit);//set it to be the end of the file and return it
			return is;
		default://if none of the above, continue on.
			t += (char)c;
		}
	}

}

int main(int argc, char* argv[])
{
	int i, count;
	double percent;
	string structure_name, out_name, line, label, perc, tmp, query,tempstring1, tempstring2;
	ifstream structure;
	ofstream out;	
	bool interactive_mode = false;

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nparse_structure_output:\n";
			cout << "Extract the cluster data from structure output\n";
			cout << "-f:\tInput Filename (include path)\n";
			cout << "-h:\tPrints this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
		else
			interactive_mode = true;
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nparse_structure_output:\n";
			cout << "Extract the cluster data from structure output\n";
			cout << "-f:\tInput Filename (include path)\n";
			cout << "-h:\tPrints this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-f")
			structure_name = tempstring2;
	}

	
	if (interactive_mode)
	{
		cout << "\nWhat is the filename?\n";
		cin >> structure_name;
		cout << "\nFilename is " << structure_name;
		cout << "\n\nProceed? (y to proceed)\n";
		cin >> query;

		if (query != "y" && query != "Y")
		{
			cout << "\n\nEnter an integer to exit!!\n";
			cin >> i;
			return 0;
		}
	}

	cout << "\nInput file name: " << structure_name << '\n';
	cout << "\n\nProceeding...\n";
	out_name = structure_name + "_clusters.txt";

	structure.open(structure_name);
	FileTest(structure, structure_name);
	out.open(out_name);
	count = 0;
	while (universal_getline(structure, line))
	{
		if (!structure.eof())
		{
			if (line == "Inferred ancestry of individuals:")
				count = 1;
			if (line == "")
				count = 0;
			if (count > 0)
			{
				if (count > 2)
				{
					stringstream ss;
					ss.str(line);
					ss >> i >> label >> perc >> tmp;
					if (count == 3)
						out << label;
					else
						out << '\n' << label;
					while (ss >> percent)
						out << '\t' << percent;
				}
				count++;
			}
		}
	}
	structure.close();
	out.close();
	if (interactive_mode)
	{
		cout << "\nDone. Input integer to quit.\n";
		cin >> i;
		return 0;
	}
	else
		return 0;
}