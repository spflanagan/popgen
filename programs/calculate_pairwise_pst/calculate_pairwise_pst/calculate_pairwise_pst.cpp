//Author: Sarah P. Flanagan
//Date: 26 May 2015
//Purpose: calculate pairwise Pst (phenotypic differentiation)
//Pst = variance between/(variance between + 2*variance within)
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

class trait
{
public:
	double value;

	trait()
	{
		value = double();
	}
};

class individual
{
public:
	string ind_id;
	vector <trait> phenotype;

	individual()
	{
		ind_id = string();
		phenotype = vector <trait>();
	}
};

class population
{
public:
	string pop_id;
	vector <individual> samples;
	int sample_num;

	population()
	{
		pop_id = string();
		samples = vector <individual>();
		sample_num = 0;
	}
};

class pst
{
public:
	string pop1, pop2;
	vector <trait> psts;

	pst()
	{
		pop1 = string();
		pop2 = string();
		psts = vector<trait>();
	}

	double calculate_pairwise_pst(population &data, population &data2, int which_trait)
	{
		//Variance components calculated using model I from Sokal and Rohlf 1981, Biometry (pp. 207-218).
		//Calculate sums of squares and then use mean squares as variance components in Pst calculation.
		//Pst calculated same as Qst (Spitze 1993), adapted from Fst (Wright).
		size_t tt;
		double grand_total, sum_sq_total, sum_sq_grp1, sum_sq_grp2;
		double correction_term, ss_groups, ss_within, ss_total;
		double between, within, pst;
		int n_total, n_groups, n_grp1, n_grp2;
		sum_sq_total = grand_total = 0;
		sum_sq_grp1 = sum_sq_grp2 = 0;
		n_total = n_grp1 = n_grp2 = 0;
		n_groups = 2;
		for (tt = 0; tt < data.samples.size(); tt++)
		{
			grand_total = grand_total + data.samples[tt].phenotype[which_trait].value;
			sum_sq_total = sum_sq_total + data.samples[tt].phenotype[which_trait].value*data.samples[tt].phenotype[which_trait].value;
			sum_sq_grp1 = sum_sq_grp1 + data.samples[tt].phenotype[which_trait].value;
			n_total++;
			n_grp1++;
		}
		sum_sq_grp1 = (sum_sq_grp1 * sum_sq_grp1)/n_grp1;
		for (tt = 0; tt < data2.samples.size(); tt++)
		{
			grand_total = grand_total + data2.samples[tt].phenotype[which_trait].value;
			sum_sq_total = sum_sq_total + data2.samples[tt].phenotype[which_trait].value*data2.samples[tt].phenotype[which_trait].value;
			sum_sq_grp2 = sum_sq_grp2 + data2.samples[tt].phenotype[which_trait].value;
			n_total++;
			n_grp2++;
		}
		sum_sq_grp2 = (sum_sq_grp2 * sum_sq_grp2)/n_grp2;
		correction_term = (grand_total * grand_total) / (n_grp1 + n_grp2);
		//cout << "Grand total:\t" << grand_total << "\nCorrection Term:\t" << correction_term;
		ss_total = sum_sq_total - correction_term;
		ss_groups = (sum_sq_grp1 + sum_sq_grp2) - correction_term;
		ss_within = ss_total - ss_groups;
		between = ss_groups / (n_groups - 1);
		within = ss_within / ((n_grp1 + n_grp2) - n_groups);
		pst = between / (between + (2*within));
		return pst;
	}
};



double calc_mean_phenotype(population &data, int phen_num)
{
	int size = 0;
	size_t t;
	double mean = 0;
	for (t = 0; t < data.samples.size(); t++)
	{
		mean = mean + data.samples[t].phenotype[phen_num].value;
		size++;
	}
	mean = mean / size;
	return mean;
}

double calculate_within_pop_variance(population &data, int phen_num)
{
	int size = 0;
	size_t t, ttt;
	double var, mean;
	var = mean = 0;
	mean = calc_mean_phenotype(data, phen_num);
	for (t = 0; t < data.samples.size(); t++)
	{
		var = var + (data.samples[t].phenotype[phen_num].value - mean)*(data.samples[t].phenotype[phen_num].value - mean);
		size++;
	}
	var = var / size;
	return var;
}

double calculate_btwn_pop_variance(population &data, population &data2, int phen_num)
{
	int size = 0;
	size_t t;
	double var, mean;
	var = mean = 0;
	mean = calc_mean_phenotype(data, phen_num);
	mean = mean + calc_mean_phenotype(data2, phen_num);
	mean = mean / 2;
	for (t = 0; t < data.samples.size(); t++)
	{
		var = var + (data.samples[t].phenotype[phen_num].value - mean)*(data.samples[t].phenotype[phen_num].value - mean);
		size++;
	}
	for (t = 0; t < data2.samples.size(); t++)
	{
		var = var + (data2.samples[t].phenotype[phen_num].value - mean)*(data2.samples[t].phenotype[phen_num].value - mean);
		size++;
	}
	var = var / size;
	return var;
}

void standardize_phenotypes(vector <population> &pop)
{
	//So I need to standardize to a given mean and variance
	size_t s, ss, sss;
	vector <double> alleles;
	double mean, std_dev;
	int num_alleles = 0;

	//calculate mean and std dev
	for (sss = 0; sss < pop.size(); sss++)
	{
		for (s = 0; s < pop[sss].samples[0].phenotype.size(); s++)
		{
			mean = calc_mean_phenotype(pop[sss], s);
			std_dev = sqrt(calculate_within_pop_variance(pop[sss], s));
			for (ss = 0; ss < pop[sss].samples.size(); ss++)
			{
				pop[sss].samples[ss].phenotype[s].value = (pop[sss].samples[ss].phenotype[s].value - mean) / std_dev;
			}
		}
	}
}

int main()
{
	int end, count, num_traits, num_pops, pop_index, i, ii,iii;
	size_t t;
	double d_temp;
	string phenotype_file_name, line, table_suffix, pst_table_name, pop_id, ind_id, temp;
	ifstream phenotype_file;
	ofstream pst_table;	
	vector <population> pops;
	vector <string> trait_names;
	vector <pst> pst_values;
	bool found, standardize, ind_id_present;

	ind_id_present = false;
	standardize = false;
	phenotype_file_name = "fem.pheno.txt";
	table_suffix = ".fem.aov.pst.txt";

	phenotype_file.open(phenotype_file_name);
	FileTest(phenotype_file, phenotype_file_name);
	count = 0;
	while (universal_getline(phenotype_file, line))
	{
		if (!phenotype_file.eof())
		{
			if (count == 0)
			{
				//it's the header
				stringstream ss(line);
				ss >> pop_id;
				if(ind_id_present)
					ss >> ind_id;
				while (ss >> temp)
					trait_names.push_back(temp);
			}
			else
			{
				//it's data
				stringstream ss(line);
				ss >> pop_id;
				found = false;
				for (t = 0; t < pops.size(); t++)
				{
					if (pops[t].pop_id == pop_id)
					{
						pop_index = t;
						found = true;
					}
				}
				if (found == false)
				{
					pop_index = pops.size();
					pops.push_back(population());
					pops[pop_index].pop_id = pop_id;
				}
				if (ind_id_present)
					ss >> ind_id;
				pops[pop_index].samples.push_back(individual());
				pops[pop_index].samples[pops[pop_index].sample_num].ind_id = ind_id; 
				for (t = 0; t < trait_names.size(); t++)
				{
					ss >> d_temp;
					pops[pop_index].samples[pops[pop_index].sample_num].phenotype.push_back(trait());
					pops[pop_index].samples[pops[pop_index].sample_num].phenotype[t].value = d_temp;
				}
				pops[pop_index].sample_num++;
			}

			count++;
		}
	}
	phenotype_file.close();
	num_pops = pops.size();
	num_traits = trait_names.size();
	cout << "\nFound " << pops.size() << " populations.\n";
	
	if (standardize)
	{
		cout << "\nStandardizing trait values.\n";
		standardize_phenotypes(pops);
	}
	int pst_index = 0;
	for (i = 0; i < num_pops - 1; i++)
	{
		for (ii = i + 1; ii < num_pops; ii++)
		{
			pst_values.push_back(pst());
			pst_values[pst_index].pop1 = pops[i].pop_id;
			pst_values[pst_index].pop2 = pops[ii].pop_id;
			for (iii = 0; iii < num_traits; iii++)
			{
				pst_values[pst_index].psts.push_back(trait());
				pst_values[pst_index].psts[iii].value = pst_values[pst_index].calculate_pairwise_pst(pops[i], pops[ii], iii);
			}
			
			pst_index++;
		}
	}

	for (iii = 0; iii < num_traits; iii++)
	{
		pst_table_name = trait_names[iii] + table_suffix;
		pst_table.open(pst_table_name);
		cout << "\nWriting to " << pst_table_name << ".\n";
		for (i = 0; i < num_pops; i++)
			pst_table << '\t' << pops[i].pop_id;
		pst_index = 0;
		for (i = 0; i < num_pops - 1; i++)
		{
			pst_table << '\n' << pops[i].pop_id << '\t';
			for (int x = 0; x < i; x++)
				pst_table << '\t';
			for (ii = i + 1; ii < num_pops; ii++)
			{
				pst_table << '\t' << pst_values[pst_index].psts[iii].value;
				pst_index++;
			}
		}
		pst_table.close();
	}

	cout << "\nDone! Input integer to quit.\n";
	cin >> end;
	return 0;
}