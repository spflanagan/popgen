//Date: 28 August 2014
//Author: Sarah P. Flanagan, sflanagan@bio.tamu.edu
//Purpose: to take two files with lists of significant scaffolds and compare them.
//this will determine how much overlap exists between population comparisons. 
//This program outputs another text file with the names of the matching scaffolds/contigs in it

//Usage:
//compare_sig_scaffolds
//compares two lists of scaffold names and outputs a file with all of those that match
//-a: first list of scaffolds (include path)
//-b: second list of scaffolds (include path)
//-o: output file name (include path)
//-h: display this message
//no arguments: interactive mode

//this program can either be run in interactive mode (responding to prompts) or by running it with the appropriate flags.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

void FileTest(ifstream& file, string filename)
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

}

int main(int argc, char* argv[])
{
	int end, count, i, ii;
	vector <string> list1, list2;
	string list1_name, list2_name, out_name, line;
	ifstream list1_file, list2_file;
	ofstream out_file;


	bool interactivemode = false;

	string query;
	string tempstring1, tempstring2, tempstring3;
	list1_name = list2_name = "default";
	out_name = "compare.txt";

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\ncompare_sig_scaffolds:\n";
			cout << "compares two lists of scaffold names and outputs a file with all of those that match\n";
			cout << "-a:\tfirst list of scaffolds (include path)\n";
			cout << "-b:\tsecond list of scaffolds (include path)\n";
			cout << "-o:\toutput file name (include path)\n";
			cout << "-h:\tdisplay this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\ncompare_sig_scaffolds:\n";
			cout << "compares two lists of scaffold names and outputs a file with all of those that match\n";
			cout << "-a:\tfirst list of scaffolds (include path)\n";
			cout << "-b:\tsecond list of scaffolds (include path)\n";
			cout << "-o:\toutput file name (include path)\n";
			cout << "-h:\tdisplay this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-a")
			list1_name = tempstring2;
		if (tempstring1 == "-b")
			list2_name = tempstring2;
		if (tempstring1 == "-o")
			out_name = tempstring3;
	}

	if (list1_name == "default")
	{
		cout << "file name of the first list of scaffolds (including path):\n";
		cin >> list1_name;
		interactivemode = true;
	}

	if (list2_name == "default")
	{
		cout << "\nfile name of the second list of scaffolds (including path):\n";
		cin >> list2_name;
		interactivemode = true;
	}

	if (out_name == "compare")
	{
		cout << "\nOutput name (include path). \nIf the default (compare.txt, in the folder in which this program is located) is acceptable, enter the character d\n";
		cin >> out_name;
		if (out_name == "d")
			out_name = "compare.txt";
		interactivemode = true;
	}

	cout << "\n\nfirst list:\t" << list1_name;
	cout << "\nsecond list:\t" << list2_name;
	cout << "\noutput name:\t" << out_name;

	if (interactivemode)
	{
		cout << "\n\nProceed? (y to proceed)\n";
		cin >> query;

		if (query != "y" && query != "Y")
		{
			cout << "\n\nEnter an integer to exit!!\n";
			cin >> i;
			return 0;
		}
	}

	cout << "\n\nProceeding...\n";


	list1_file.open(list1_name);
	FileTest(list1_file, list1_name);
	
	while (getline(list1_file, line))
	{
		list1.push_back(line);
	}
	list1_file.close();

	list2_file.open(list2_name);
	FileTest(list2_file, list2_name);
	
	while (getline(list2_file, line))
	{
		list2.push_back(line);
	}
	list1_file.close();

	out_file.open(out_name);
	count = 0;
	for (i = 0; i < list1.size(); i++)
	{
		for (ii = 0; ii < list2.size(); ii++)
		{
			if (list1[i] == list2[ii])
			{
				count++;
				out_file << list1[i] << '\n';
			}
		}
	}
	out_file.close();
	cout << "Found " << count << " matching scaffolds!\n";

	if (interactivemode)
	{
		cout << "Done! Enter integer to quit.\n";
		cin >> end;
	}

	return 0;
}