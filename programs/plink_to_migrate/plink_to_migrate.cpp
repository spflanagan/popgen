//Author: Sarah P. Flanagan
//Date: 9 March 2016
//Purpose: To convert SNP data in PLINK format into MIGRATE-N format
//Input: Plink map and ped files
//Output: A migrate file in HapMap format.
//The "Family ID" field in the Plink ped file should be the population name

//Plink format: two files, a map file and a ped file
//map file has four columns: Chr, SNPID, MapDistance, BP
//ped file has Family ID,Individual ID,Paternal ID,Maternal ID,Sex(1 = male; 2 = female; other = unknown),Phenotype (-9 is missing), 
//followed by two columns for each locus.

//migrate format has:
//<Number of populations> <number of loci>[project title 0 - 79]
//<Any Number> <title for population 0 - 79>
//<Position on chromosome locus1> <TAB><allele><TAB><number><TAB><allele><TAB><number><TAB><total>
//<Position on chromosome locus2> <TAB><allele><TAB><number><TAB><allele><TAB><number><TAB><total>
//....
//<Any Number> <title for population 0 - 79>
//<Position on chromosome locus1> <TAB><allele><TAB><number><TAB><allele><TAB><number><TAB><total>
//<Position on chromosome locus2> <TAB><allele><TAB><number><TAB><allele><TAB><number><TAB><total>

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

class reference
{
public:
	string id;
	vector<string> alleles;

	reference()
	{
		id = string();
		alleles = vector<string>();
	}

};

class locus
{
public:
	vector<string> alleles;
	vector<int> allele_counts;

	locus()
	{
		alleles = vector<string>();
		allele_counts = vector<int>();
	}
};

class population
{
public:
	int pop_num, ind_count;
	string pop_name;
	vector<locus> loci;

	population()
	{
		pop_num = int();
		ind_count = int();
		pop_name = string();
		loci = vector<locus>();
	}
};

int main()
{
	int i, ii, iii, num_pops, num_loci, index, count, all1i, all2i;
	string chr, pos, id, dist;
	string famid, indid, patid, matid, sex, phen, allele1, allele2;
	string ped_name, map_name, migrate_name,line;
	vector<reference> loci;
	
	ifstream ped, map;
	ofstream migrate;
	vector<population> pops;

	ped_name = "../../sw_results/migrate/subset.ped"; //this has been converted to have the pop IDs in Fam ID column.
	map_name = "../../sw_results/stacks/populations/subset.map";
	migrate_name = "../../sw_results/migrate/migrate_in.txt";

	//read in the map file first
	map.open(map_name);
	FileTest(map, map_name);
	num_loci = 0;	
	while (universal_getline(map, line))
	{
		if (!map.eof())
		{
			if (line.substr(0, 1) != "#")
			{
				stringstream ss;
				ss.str(line);
				ss >> chr >> id >> dist >> pos;
				loci.push_back(reference());
				loci.back().id = id;
				num_loci++;
			}
		}
	}
	map.close();

	//read in the ped file
	ped.open(ped_name);
	FileTest(ped, ped_name);
	num_pops = 0;
	while (universal_getline(ped, line))
	{
		if (!ped.eof())
		{
			if (line.substr(0, 1) != "#")
			{
				stringstream ss;
				ss.str(line);
				ss >> famid >> indid >> patid >>  matid >> sex >> phen;
				index = -5;
				for (i = 0; i < pops.size(); i++)
				{
					if (pops[i].pop_name == famid)
						index = i;
				}
				if (index < 0)
				{
					pops.push_back(population());
					index = pops.size() - 1;
					pops[index].pop_name = famid;
					pops[index].pop_num = index;
					pops[index].ind_count = 0;
					num_pops++;
				}
				pops[index].ind_count++;
				count = 0;
				if (pops[index].loci.size() == 0)//then it's the first individual in the population
				{
					while (ss >> allele1 >> allele2)
					{
						pops[index].loci.push_back(locus());
						if (allele1 != "0")
						{
							//see if they're in the reference already.
							all1i = all2i = -5;
							for (i = 0; i < loci[count].alleles.size(); i++)
							{
								if (loci[count].alleles[i] == allele1)
									all1i = i;
								if (loci[count].alleles[i] == allele2)
									all2i = i;

							}
							if (all1i < 0)
							{
								loci[count].alleles.push_back(allele1);
								all1i = loci[count].alleles.size() - 1;
							}

							if (all2i < 0 && allele1 != allele2)
							{
								loci[count].alleles.push_back(allele2);
								all2i = loci[count].alleles.size() - 1;
							}

							//now add it to the population counts
							if (allele1 == allele2)
							{
								pops[index].loci[count].alleles.push_back(loci[count].alleles[all1i]);
								pops[index].loci[count].allele_counts.push_back(2);
							}
							else
							{
								pops[index].loci[count].alleles.push_back(loci[count].alleles[all1i]);
								pops[index].loci[count].allele_counts.push_back(1);
								pops[index].loci[count].alleles.push_back(loci[count].alleles[all2i]);
								pops[index].loci[count].allele_counts.push_back(1);
							}
						}
						count++;
					}
				}
				else
				{
					while (ss >> allele1 >> allele2)
					{//find them in the index
						if (allele1 != "0")
						{
							//see if they're in the reference already.
							all1i = all2i = -5;
							for (i = 0; i < loci[count].alleles.size(); i++)
							{
								if (loci[count].alleles[i] == allele1)
									all1i = i;
								if (loci[count].alleles[i] == allele2)
									all2i = i;

							}
							if (all1i < 0)
							{
								loci[count].alleles.push_back(allele1);
								all1i = loci[count].alleles.size() - 1;
							}

							if (all2i < 0 && allele1 != allele2)
							{
								loci[count].alleles.push_back(allele2);
								all2i = loci[count].alleles.size() - 1;
							}
							
							//now check the pop
							bool found1, found2;
							found1 = found2 = false;
							for (i = 0; i < pops[index].loci[count].alleles.size(); i++)
							{
								if (pops[index].loci[count].alleles[i] == allele1)
								{
									pops[index].loci[count].allele_counts[i]++;
									found1 = true;
								}
								if (pops[index].loci[count].alleles[i] == allele2)
								{
									pops[index].loci[count].allele_counts[i]++;
									found2 = true;
								}
							}
							if (!found1)
							{
								pops[index].loci[count].alleles.push_back(allele1);
								pops[index].loci[count].allele_counts.push_back(1);
							}
							if (!found2 && allele1 != allele2)
							{
								pops[index].loci[count].alleles.push_back(allele2);
								pops[index].loci[count].allele_counts.push_back(1);
							}
						}
						count++;
					}
				}//else it's not the first ind in the pop
			}
		}
	}//while ped
	ped.close();
	//make sure both snps are represented
	for (i = 0; i < pops.size(); i++)
	{
		for (ii = 0; ii < loci.size(); ii++)
		{
			if (pops[i].loci[ii].alleles.size() != loci[ii].alleles.size())
			{
				vector<bool> found;
				for (iii = 0; iii < loci[ii].alleles.size(); iii++)
				{
					found.push_back(false);
					for (int j = 0; j < pops[i].loci[ii].alleles.size(); j++)
					{
						if (loci[ii].alleles[iii] == pops[i].loci[ii].alleles[j])
							found[iii] = true;							
					}
				}
				for (iii = 0; iii < loci[ii].alleles.size(); iii++)
				{
					if (!found[iii])
					{
						pops[i].loci[ii].alleles.push_back(loci[ii].alleles[iii]);
						pops[i].loci[ii].allele_counts.push_back(0);
					}
				}
			}
		}
	}
	cout << "\nParsed " << num_loci << " loci in " << num_pops << " populations.\n";

	cout << "\nWriting loci to Migrate file\n";
	migrate.open(migrate_name);
	migrate << num_pops << " " << num_loci;
	for (i = 0; i < num_pops; i++)
	{
		migrate << '\n' << pops[i].ind_count*2 << " POP" << pops[i].pop_name;
		for (ii = 0; ii < num_loci; ii++)
		{
			migrate << '\n' << ii;
			int sum = 0;
			for (iii = 0; iii < pops[i].loci[ii].alleles.size(); iii++)
			{
				migrate << '\t' << pops[i].loci[ii].alleles[iii] << '\t' << pops[i].loci[ii].allele_counts[iii];
				sum = sum + pops[i].loci[ii].allele_counts[iii];
			}
			migrate << '\t' << sum;
		}
	}
	migrate.close();

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}