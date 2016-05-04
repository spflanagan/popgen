//Author: Sarah P. Flanagan
//Date: 13 June 2015
//Purpose: parse a World Oceans Database .csv output file
//extract the desired data and output in a more useful spreadsheet format

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

class wod_record
{
public:
	double latitude;
	double longitude;
	int year, month, day;
	double depth;
	double temp, salinity, oxygen, pH;

	wod_record()
	{
		latitude = longitude = double();
		year = month = day = int();
		depth = double();
		temp = salinity = oxygen = pH = -5;
	}
};

int main(int argc, char* argv[])
{
	int i, end, count, num_var;
	string csv_file_name, out_file_name, line;
	string tempstring1, tempstring2, query;
	ifstream csv_file;
	ofstream out_file;
	vector <wod_record> records;
	vector <string> var_names;
	bool interactive_mode = false;

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nparse_wod_data:\n";
			cout << "Read in a World Oceans Database .csv file and extract relevant info.\n";
			cout << "-f:\tFilename (include path)\n";
			cout << "-h:\tPrints this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nparse_wod_data:\n";
			cout << "Read in a World Oceans Database .csv file and extract relevant info.\n";
			cout << "-f:\tFilename (include path)\n";
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
			csv_file_name = tempstring2;
	}

	cout << "\nInput file name: " << csv_file_name << '\n';
	if (interactive_mode)
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

	out_file_name = csv_file_name.substr(0,csv_file_name.length() - 4) + "_out.txt";
	csv_file.open(csv_file_name);
	FileTest(csv_file, csv_file_name);
	while (universal_getline(csv_file, line))
	{
		if (!csv_file.eof())
		{
			if (line.substr(0, 1) == "#")
				records.push_back(wod_record());
			else
			{
				stringstream ss(line);
				string line_comp;
				vector <string> temp;
				while (getline(ss, line_comp, ','))
					temp.push_back(line_comp);
				if (temp[0].substr(0,8) == "Latitude")
					records.back().latitude = atof(temp[2].c_str());
				if (temp[0].substr(0,9) == "Longitude")
					records.back().longitude = atof(temp[2].c_str());
				if (temp[0].substr(0,4) == "Year")
					records.back().year = atoi(temp[2].c_str());
				if (temp[0].substr(0,5) == "Month")
					records.back().month = atoi(temp[2].c_str());
				if (temp[0].substr(0,3) == "Day")
					records.back().day = atoi(temp[2].c_str());
				if (temp[0].substr(0,8) == "VARIABLE")
				{
					num_var = (temp.size() - 1) / 3;
					for (size_t t = 1; t < temp.size() - 1; t += 3)
						var_names.push_back(temp[t]);					
				}
				if (temp[0].substr(temp[0].size() -2, temp[0].size()) == " 1")
				{//it's the first reading!
					count = 0;
					for (size_t t = 1; t < temp.size() - 1; t += 3)
					{
						if (var_names[count].substr(0,5) == "Depth")
							records.back().depth = atof(temp[t].c_str());
						if (var_names[count].substr(0,4) == "Temp")
							records.back().temp = atof(temp[t].c_str());
						if (var_names[count].substr(0,8) == "Salinity")
							records.back().salinity = atof(temp[t].c_str());
						if (var_names[count].substr(0,7) == "Oxygen")
							records.back().oxygen = atof(temp[t].c_str());
						if (var_names[count].substr(0,2) == "pH")
							records.back().pH = atof(temp[t].c_str());
						count++;
					}
					

				}
			}
		}
	}
	csv_file.close();

	cout << "\nwriting to output file " << out_file_name << "\n";
	out_file.open(out_file_name);
	out_file << "Lat\tLong\tYear\tMonth\tDay\tDepth\tTemp\tSalinity\tOxygen\tpH";
	for (size_t t = 0; t < records.size(); t++)
	{
		out_file << '\n' << records[t].latitude << '\t' << records[t].longitude << '\t'
			<< records[t].year << '\t' << records[t].month << '\t' << records[t].day;
		if (records[t].depth != -5)
			out_file << '\t' << records[t].depth;
		else
			out_file << '\t' << "NA";
		if(records[t].temp != -5)
			out_file << '\t' << records[t].temp;
		else
			out_file << '\t' << "NA";
		if (records[t].salinity != -5)
			out_file << '\t' << records[t].salinity;
		else
			out_file << '\t' << "NA";
		if (records[t].oxygen != -5)
			out_file << '\t' << records[t].oxygen;
		else
			out_file << '\t' << "NA";
		if (records[t].pH != -5)
			out_file << '\t' << records[t].pH;
		else
			out_file << '\t' << "NA";
	}
	out_file.close();

	if (interactive_mode)
	{
		cout << "\nDone! Input integer to quit.\n";
		cin >> end;
		return 0;
	}
	else
	{
		cout << "\nDone!\n";
		return 0;
	}
}