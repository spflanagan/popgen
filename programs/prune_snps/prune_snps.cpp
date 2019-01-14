//Author: Sarah P. Flanagan
//Date: 3 March 2016
//Purpose: To take a plink map file and prune it for one snp per locus and one locus per 2kb. Then prune a ped file to match the updated map file. 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "../random_numbers.h"

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

string find_and_replace(string &s, string toReplace, string replaceWith)
{
	size_t pos = 0;
	while ((pos = s.find(toReplace, pos)) != string::npos) {
		s.replace(pos, toReplace.length(), replaceWith);
		pos += replaceWith.length();
	}
	return s;
}

class locus_information
{
public:
	vector<int> bp, indices;
	string chromosome, locus_id;
	vector<string> snp_id;
	
	locus_information()
	{
		bp = indices = vector<int>();
		chromosome = locus_id = string();
		snp_id = vector <string>();
	}
};
int main(int argc, char* argv[])
{
	int i, ii,iii, m,bp, index, line_count,count,loc_id; 
	vector<int> keep_index; //give it a 0 if it's removed but a 1 if it's kept.
	string chr, snp, locid;
	string fam_id, ind_id, pat_id, mat_id, sex, phenotype, all1, all2;
	string line, map_name, ped_name, out_map_name, out_ped_name;
	ifstream map, ped;
	ofstream out_map, out_ped;
	vector<locus_information> loci;
	string query;
	string tempstring1, tempstring2, tempstring3;
	bool interactivemode;
	ped_name = "plink.ped";
	map_name = "plink.map";
	out_ped_name = "pruned.ped";
	out_map_name = "pruned.map";

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nprune_snps:\n";
			cout << "Prunes PLINK files to remove multiple SNPs per RAD locus\nAlso prunes RAD loci too close to each other.\n";
			cout << "-p:\tInput Plink .ped file (include path)\n";
			cout << "-m:\tInput plink .map file (include path)\n";
			cout << "-o:\toutput name (which will be used for both ped and map files)\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
		else
			interactivemode = true;
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nprune_snps:\n";
			cout << "Prunes PLINK files to remove multiple SNPs per RAD locus\nAlso prunes RAD loci too close to each other.\n";
			cout << "-p:\tInput Plink .ped file (include path)\n";
			cout << "-m:\tInput plink .map file (include path)\n";
			cout << "-o:\toutput name (which will be used for both ped and map files)\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-p")
			ped_name = tempstring2;
		if (tempstring1 == "-m")
			map_name = tempstring2;
		if (tempstring1 == "-o")
			out_ped_name = tempstring2;
	}
	out_map_name = out_ped_name + ".map";
	out_ped_name = out_ped_name + ".ped";

	if (interactivemode)
	{
		cout << "\nProvide input ped file name and path.\n";
		cin >> ped_name;
		cout << "\nProvide input map file name and path.\n";
		cin >> map_name;
		cout << "\nProvide output ped file name and path.\n";
		cin >> out_ped_name;
		cout << "\nProvide output map file name and path.\n";
		cin >> out_map_name;
		cout << "\nInput names: " << ped_name << ", " << map_name;
		cout << "\nOutput names: " << out_ped_name << ", " << out_map_name;
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

	map.open(map_name);
	FileTest(map, map_name);
	count = line_count = 0; 
	while (universal_getline(map, line))
	{
		if (!map.eof())
		{
			if (line.substr(0, 1) != "#")
			{
				stringstream ss;
				ss.str(line);
				ss >> chr >> snp >> m >> bp; //plink map file has Chromosome\tSNPID\tdistance in Morgans\tBP
				stringstream ssnp;
				ssnp << snp;
				getline(ssnp, locid, '_');
				index = -5;
				for (i = 0; i < loci.size(); i++)
				{
					if (loci[i].locus_id == locid)
						index = i;
				}
				if (index >= 0)
				{
					loci[index].bp.push_back(bp);
					loci[index].snp_id.push_back(snp);
					loci[index].indices.push_back(count);
				}
				else
				{
					loci.push_back(locus_information());
					loci.back().chromosome = chr;
					loci.back().locus_id = locid;
					loci.back().bp.push_back(bp);
					loci.back().snp_id.push_back(snp);
					loci.back().indices.push_back(count);
				}
				keep_index.push_back(1);
				count++;
			}
			line_count++;
			if (line_count % 1000 == 0)
				cout << "\nProcessed " << line_count << " SNP identifiers.";
		}
	}
	map.close();
	cout << "\nThis dataset contains " << line_count << " SNPs.";

	cout << "\nNow keeping only one SNP per RAD locus.";
	for (i = 0; i < loci.size(); i++)
	{
		if (loci[i].snp_id.size() > 1)
		{
			index = randnum(loci[i].snp_id.size());
			for (ii = 0; ii < loci[i].snp_id.size(); ii++)
			{
				if (index != ii)
					keep_index[loci[i].indices[ii]] = 0;
			}
		}
	}
	count = 0;
	for (i = 0; i < keep_index.size(); i++)
	{
		if (keep_index[i] == 1)
			count++;
	}
	cout << '\n' << count << " SNPs were retained.";

	cout << "\nRemoving loci within 2kb of each other.";
	for (i = 0; i < loci.size(); i++)
	{
		index = -1;
		for (ii = 0; ii < loci[i].indices.size() - 1; ii++)
		{
			if (keep_index[loci[i].indices[ii]] == 1)
				index = ii;
		}
		if (index >= 0)
		{
			for (ii = i+1; ii < loci.size(); ii++)
			{
				if (loci[i].chromosome == loci[ii].chromosome) //if they're on the same chromosome
				{
					loc_id = -1;
					for (iii = 0; iii < loci[ii].indices.size(); iii++) //check to see if any of the snps are kept
					{
						if (keep_index[loci[ii].indices[iii]] == 1)
							loc_id = iii;
					}
					if (loc_id >= 0 && abs(loci[i].bp[index] - loci[ii].bp[loc_id]) < 1000)
					{
						if (genrand() < 0.5)
							keep_index[loci[i].indices[index]] = 0;
						else
							keep_index[loci[ii].indices[loc_id]] = 0;
					}
				}
			}
		}//if index > - 0
	}

	cout << "\nRe-writing map file.";
	out_map.open(out_map_name);
	count = 0;
	for (i = 0; i < loci.size(); i++)
	{
		for (ii = 0; ii < loci[i].indices.size(); ii++)
		{
			if (keep_index[loci[i].indices[ii]] == 1)
			{
				if (count == 0)
					out_map << loci[i].chromosome << '\t' << loci[i].snp_id[ii] << "\t0\t" << loci[i].bp[ii];
				else
					out_map << '\n' << loci[i].chromosome << '\t' << loci[i].snp_id[ii] << "\t0\t" << loci[i].bp[ii];
				count++;
			}
		}
	}
	out_map.close();
	cout << '\n' << count << " SNPs were retained.";

	ped.open(ped_name);
	FileTest(ped, ped_name);
	out_ped.open(out_ped_name);
	cout << "\nRe-writing ped file to " << out_ped_name;
	line_count = 0;
	while (universal_getline(ped, line))
	{
		if (!ped.eof())
		{
			if (line.substr(0, 1) != "#")
			{
				stringstream ss;
				ss.str(line);
				ss >> fam_id >> ind_id >> pat_id >> mat_id >> sex >> phenotype;
				if (line_count == 0)
					out_ped << fam_id << '\t' << ind_id << '\t' << pat_id << '\t' << mat_id << '\t' << sex << '\t' << phenotype;
				else
					out_ped << '\n' << fam_id << '\t' << ind_id << '\t' << pat_id << '\t' << mat_id << '\t' << sex << '\t' << phenotype;
				count = 0;
				while (ss >> all1 >> all2)
				{
					if (keep_index[count] == 1)
						out_ped << '\t' << all1 << '\t' << all2;
					count++;
				}
				line_count++;
			}
		}
	}
	ped.close();
	out_ped.close();

	cout << "\nDone! Enter integer to quit.\n";
	cin >> i;
	return 0;
}