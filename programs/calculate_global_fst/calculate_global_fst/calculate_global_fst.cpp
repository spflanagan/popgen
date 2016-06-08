//Author: Sarah P. Flanagan
//Date: 9 June 2015
//Purpose: calculate global Fst 
//Using Nei's formulation:
//Fst=1-((sum of each population's expected heterozygosity)/((number of populations)*(overall expected heterozygosity)))
//input = plink files, no header on ped file, and each locus has two columns per ind
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>

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

class locus{
public:
	string id;
	int bp;
	string scaffold;
	vector <double> allele_freq;
	double fst;
	double he;
	double ho;

	locus()
	{
		id = string();
		bp = int();
		scaffold = string();
		allele_freq = { 0, 0, 0, 0 };
		fst = double();
		he = 1;
		ho = 0;
	}

	int determine_allele_locus(char nt)
	{
		int nt_loc = 5;
		vector <char> ch_nt = { 'A', 'C', 'T', 'G' };
		for (int i = 0; i < 4; i++)
		{
			if (ch_nt[i] == nt)
				nt_loc = i;
		}
		return(nt_loc);
	}
};


class population
{
public:
	string pop_id;
	vector <locus> genotypes;
	int num_inds;

	population()
	{
		pop_id = string();
		genotypes = vector <locus>();
		num_inds = 0;
	}

};


int main()
{
	int i, ii, end, count, num_pops, pop_index, loc_id, mat, pat, sex, phen, dist, bp, gen_index;
	char gen1, gen2;
	size_t t, tt, ttt;
	string pop, ind, chrom, snp, line, pop_id, temp;
	string out_file_name, ped_name, map_name;
	ifstream ped_file, map_file;
	ofstream out_file;
	vector <population> pops;
	population overall;
	bool found;

	
	ped_name = "B://ubuntushare//popgen//nerophis//stacks//pruned.ped";
	map_name = "B://ubuntushare//popgen//nerophis//stacks//pruned.map";
	out_file_name = "B://ubuntushare//popgen//nerophis//stacks//pruned.globalstats.txt";

	ped_file.open(ped_name);
	FileTest(ped_file, ped_name);
	count = 0;
	while (universal_getline(ped_file, line))
	{
		if (!ped_file.eof())
		{
			if (count == 0)
			{

				stringstream ss(line);
				ss >> pop >> ind >> mat >> pat >> sex >> phen;
				pops.push_back(population());
				pops.back().pop_id = pop;
				pop_index = pops.size() - 1;
				while (ss >> gen1 >> gen2)
				{
					pops[pop_index].genotypes.push_back(locus());
					overall.genotypes.push_back(locus());
					if (gen1 != '0' && gen2 != '0')
					{
						pops[pop_index].genotypes.back().allele_freq[pops[pop_index].genotypes.back().determine_allele_locus(gen1)]++;
						pops[pop_index].genotypes.back().allele_freq[pops[pop_index].genotypes.back().determine_allele_locus(gen2)]++;
						overall.genotypes.back().allele_freq[overall.genotypes.back().determine_allele_locus(gen1)]++;
						overall.genotypes.back().allele_freq[overall.genotypes.back().determine_allele_locus(gen2)]++;
						if (gen1 == gen2)
						{
							pops[pop_index].genotypes.back().ho++;
							overall.genotypes.back().ho++;
						}
					}
				}
				count++;
			}
			else
			{
				stringstream ss(line);
				ss >> pop >> ind >> mat >> pat >> sex >> phen;
				found = false;
				for (t = 0; t < pops.size(); t++)
				{
					if (pops[t].pop_id == pop)
					{
						pop_index = t;
						found = true;
					}
				}
				if (found == false)
				{
					pops.push_back(population());
					pops.back().pop_id = pop;
					pop_index = pops.size() - 1;
					for (t = 0; t < pops[0].genotypes.size(); t++)
						pops[pop_index].genotypes.push_back(locus());
				}
				
				pops[pop_index].num_inds++;
				overall.num_inds++;
				gen_index = 0;
				while (ss >> gen1 >> gen2)
				{//0=A,1=C,2=T,3=G
					if (gen1 != '0' && gen2 != '0')
					{
						pops[pop_index].genotypes[gen_index].allele_freq[pops[pop_index].genotypes.back().determine_allele_locus(gen1)]++;
						pops[pop_index].genotypes[gen_index].allele_freq[pops[pop_index].genotypes.back().determine_allele_locus(gen2)]++;
						overall.genotypes[gen_index].allele_freq[overall.genotypes.back().determine_allele_locus(gen1)]++;
						overall.genotypes[gen_index].allele_freq[overall.genotypes.back().determine_allele_locus(gen2)]++;
						if (gen1 == gen2)
						{
							pops[pop_index].genotypes[gen_index].ho++;
							overall.genotypes[gen_index].ho++;
						}
					}
					gen_index++;
				}
				count++;
			}
		}
	}
	ped_file.close();
	cout << "\nFinished parsing ped file.\n";

	
	map_file.open(map_name);
	FileTest(map_file, map_name);
	count = 0;
	while (universal_getline(map_file, line))
	{
		if (!map_file.eof())
		{
			stringstream ss(line);
			ss >> chrom >> snp >> dist >> bp;
			overall.genotypes[count].bp = bp;
			overall.genotypes[count].scaffold = chrom;
			overall.genotypes[count].id = snp;
			count++;
			
		}
	}
	map_file.close();
	cout << "\nFinished parsing map file\n";

	//need to get the allele *frequencies* and calculate exp het (1-sum(f*f)), where f = freq of each allele
	for (t = 0; t < overall.genotypes.size(); t++)
	{
		overall.genotypes[t].he = 1;
		overall.genotypes[t].ho = overall.genotypes[t].ho / (overall.num_inds * 2);
		for (tt = 0; tt < pops.size(); tt++)
		{
			pops[tt].genotypes[t].ho = pops[tt].genotypes[t].ho / (pops[tt].num_inds * 2);
			pops[tt].genotypes[t].he = 1;
		}
		for (i = 0; i < 4; i++)
		{
			overall.genotypes[t].allele_freq[i] = overall.genotypes[t].allele_freq[i] / (overall.num_inds * 2);
			overall.genotypes[t].he = overall.genotypes[t].he - (overall.genotypes[t].allele_freq[i] * overall.genotypes[t].allele_freq[i]);
			for (tt = 0; tt < pops.size(); tt++)
			{
				pops[tt].genotypes[t].allele_freq[i] = pops[tt].genotypes[t].allele_freq[i] / (pops[tt].num_inds * 2);
				pops[tt].genotypes[t].he = pops[tt].genotypes[t].he - (pops[tt].genotypes[t].allele_freq[i] * pops[tt].genotypes[t].allele_freq[i]);
			}
		}
		overall.genotypes[t].fst = 0;
		for (tt = 0; tt < pops.size(); tt++)
		{
			overall.genotypes[t].fst = overall.genotypes[t].fst + pops[tt].genotypes[t].he;
			cout << overall.genotypes[t].fst << '\t';
		}
		overall.genotypes[t].fst = 1 - (overall.genotypes[t].fst / (pops.size() * overall.genotypes[t].he));
	}

	cout << "\nWriting Summary Statistics to file.\n";
	out_file.open(out_file_name);
	out_file << "Chrom\tLocus\tBP\tObs.Het\tExp.Het\tFst";
	for (t = 0; t < overall.genotypes.size(); t++)
	{
		out_file << '\n';
		out_file << overall.genotypes[t].scaffold << '\t' << overall.genotypes[t].id << '\t' << overall.genotypes[t].bp << '\t';
		out_file << overall.genotypes[t].ho << '\t' << overall.genotypes[t].he << '\t' << overall.genotypes[t].fst;
	}
	out_file.close();
	cout << "\nDone! Input integer to quit.\n";
	cin >> end;
	return 0;
}