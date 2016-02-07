//Date: 15 August 2014
//Author: Sarah P. Flanagan
//Purpose: take output from stacks and calculate coverage data per individual and per locus
//requires: (1) the pathname to the stacks files;
//and (2) a text file in the same directory as this program with a list of the root filenames (before .matches.tsv).

//Usage:
//coverage_calcs
//calculates coverage summary stats from individuals' matches.tsv files, outpu by stacks
//-p: path to directory with stacks *.matches.tsv file
//-i: List of sample IDs (the component of the filenames before the .matches.tsv)
//-n: output name (e.g. pop ID, without any extensions; default = coverage)
//no arguments: interactive mode

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;


class summary
{
public:
	string name;
	int stack_depth;
	int num_with_this;//number of loci for individuals, number of individuals for loci
	int num_above_5x;
	int num_above_10x;
	double average_cov;
	int min_cov;
	int max_cov;

	//Constructor definition
	summary()
	{
		name = "sample";
		stack_depth = 0;
		num_with_this = 0;
		num_above_10x = 0;
		num_above_5x = 0;
		average_cov = 0;
		min_cov = 10;
		max_cov = 0;
	}
};

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
	vector <summary> per_ind;
	vector <summary> per_loc;
	vector <summary> overall;
	string path_name, line, full_name;
	string samplelist_filename;
	ifstream samplelist, matches_file;
	string samplename, matches_name;
	int count = 0;
	int end, i, a, b, c, d, e, g;
	double h;
	string f;
	istringstream lin;
	ofstream ind_sum, loc_sum, all_sum;
	string out_name = "coverage";
	bool interactivemode = false;

	string query;
	string tempstring1, tempstring2, tempstring3;
	path_name = samplelist_filename = "default";

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\ncoverage_calcs:\n";
			cout << "calculates coverage summary stats from individuals' matches.tsv files, outpu by stacks\n";
			cout << "-p:\tpath to directory with stacks *.matches.tsv file\n";
			cout << "-i:\tList of sample IDs (the component of the filenames before the .matches.tsv)\n";
			cout << "-n:\toutput name (e.g. pop ID, without any extensions; default = coverage)\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\ncoverage_calcs:\n";
			cout << "calculates coverage summary stats from individuals' matches.tsv files, outpu by stacks\n";
			cout << "-p:\tpath to directory with stacks *.matches.tsv file\n";
			cout << "-i:\tList of sample IDs (the component of the filenames before the .matches.tsv)\n";
			cout << "-n:\toutput name (e.g. pop ID, without any extensions; default = coverage)\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-p")
			path_name = tempstring2;
		if (tempstring1 == "-i")
			samplelist_filename = tempstring2;
		if (tempstring1 == "-n")
			out_name = tempstring3;
	}

	if (path_name == "default")
	{
		cout << "path to directory with stacks *.matches.tsv file:\n";
		cin >> path_name;
		interactivemode = true;
	}

	if (samplelist_filename == "default")
	{
		cout << "\nList of sample IDs (the component of the filenames before the .matches.tsv):\n";
		cin >> samplelist_filename;
		interactivemode = true;
	}

	if (out_name == "coverage")
	{
		cout << "\nOutput name (e.g. pop ID, without any extensions). \nIf the default (coverage) is acceptable, enter the character d\n";
		cin >> out_name;
		if (out_name == "d")
			out_name = "coverage";
		interactivemode = true;
	}

	cout << "\n\npath to directory with stacks *.matches.tsv file:\t" << path_name;
	cout << "\nList of sample IDs (the component of the filenames before the .matches.tsv):\t" << samplelist_filename;

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

	//cout << "Enter path to directory with stacks *.matches.tsv files\n";
	//getline(cin, path_name, '\n');
	//path_name = "C:\Users\Sarah\Dropbox\popgen_RAD\";

	//cout << "Name of the file containing a list of the sample IDs\n";
	//getline(cin, samplelist_filename, '\n');

	samplelist.open(samplelist_filename);
	FileTest(samplelist, samplelist_filename);

	while (getline(samplelist, samplename, '\n'))
	{
		istringstream ss(samplename);
		per_ind.push_back(summary());
		ss >> per_ind[count].name;
		cout << per_ind[count].name << '\n';
		count++;
	}
	//per_loc.push_back(summary());//give it the first value
	cout << per_ind.size();
	for (i = 0; i < per_ind.size(); i++)
	{
		cout << "processing sample " << i + 1 << "\n";
		matches_name.clear();
		matches_name = path_name;
		matches_name = matches_name.append(per_ind[i].name.append(".matches.tsv"));
		//cout << matches_name << '\n';
		matches_file.open(matches_name);
		FileTest(matches_file, matches_name);
		while (getline(matches_file, line, '\n'))
		{
			lin.clear();
			lin.str(line);
			if (lin >> a >> b >> c >> d >> e >> f >> g >> h)
			{
				//sqlID,batchID,catalogID,indID,indlocusID, haplotype,stack depth, log likelihood
				per_ind[i].stack_depth = per_ind[i].stack_depth + g;
				per_ind[i].num_with_this++;
				//need to designate new locus components
				if (c >= per_loc.size())
				{
					for (int s = per_loc.size(); s < c + 1; s++)
					{
						per_loc.push_back(summary());
					}
				}
				if (c % 10000 == 0)
				{
					cout << "Working...found locus " << c << '\n';
				}
				per_loc[c].stack_depth = per_loc[c].stack_depth + g;
				per_loc[c].num_with_this++;
				per_loc[c].name = c;
				if (g >= 5)
				{
					per_loc[c].num_above_5x++;
					per_ind[i].num_above_5x++;
				}
				if (g >= 10)
				{
					per_loc[c].num_above_10x++;
					per_ind[i].num_above_10x++;
				}
				/*if (i == 0)
				{
				per_loc[c].min_cov = g;
				per_ind[i].min_cov = g;
				}*/
				if (g < per_loc[c].min_cov)
					per_loc[c].min_cov = g;
				if (g < per_ind[i].min_cov)
					per_ind[i].min_cov = g;
				if (g > per_loc[c].max_cov)
					per_loc[c].max_cov = g;
				if (g > per_ind[i].max_cov)
					per_ind[i].max_cov = g;
			}
			else
			{
				cout << "WARNING: failed to decode line '" << line << "'\n";
			}
		}
		matches_file.close();
	}//end count

	cout << "Calculating summary stats\n";
	int num_loci, num_above5x, num_above10x, overall_min, overall_max;
	double avg_per_locus, avg_per_ind;
	num_loci = num_above5x = num_above10x = overall_max = avg_per_ind = avg_per_locus = 0;
	overall_min = 1;
	full_name = out_name + "_ind_sum.txt";
	ind_sum.open(full_name.c_str());
	ind_sum << "ID\tNumLoci\tAvgCov\tMinCov\tMaxCov\t5x\t10x\tTotalStacks\n";
	for (i = 0; i < per_ind.size(); i++)
	{
		if (per_ind[i].num_with_this > 0)
			per_ind[i].average_cov = per_ind[i].stack_depth / per_ind[i].num_with_this;
		else
			per_ind[i].average_cov = 0;
		avg_per_ind = avg_per_ind + per_ind[i].average_cov;
		ind_sum << per_ind[i].name << '\t' << per_ind[i].num_with_this << '\t' << per_ind[i].average_cov << '\t' << per_ind[i].min_cov << '\t' <<
			per_ind[i].max_cov << '\t' << per_ind[i].num_above_5x << '\t' << per_ind[i].num_above_10x << '\t' << per_ind[i].stack_depth << '\n';
	}
	avg_per_ind = avg_per_ind / per_ind.size();
	ind_sum.close();

	full_name = out_name + "_loc_sum.txt";
	loc_sum.open(out_name.c_str());
	loc_sum << "ID\tNumInd\tAvgCov\tMinCov\tMaxCov\t5x\t10x\tTotalStacks\n";
	for (c = 0; c < per_loc.size(); c++)
	{
		if (per_loc[c].num_with_this > 0)
		{
			per_loc[c].average_cov = per_loc[c].stack_depth / per_loc[c].num_with_this;
			num_loci++;
		}
		else
			per_loc[c].average_cov = 0;
		if (per_loc[c].num_with_this < overall_min)
			overall_min = per_loc[c].num_with_this;
		if (per_loc[c].num_with_this > overall_max)
			overall_max = per_loc[c].num_with_this;
		if (per_loc[c].average_cov >= 5)
			num_above5x++;
		if (per_loc[c].average_cov >= 10)
			num_above10x++;
		avg_per_locus = avg_per_locus + per_loc[c].average_cov;
		loc_sum << per_loc[c].name << '\t' << per_loc[c].num_with_this << '\t' << per_loc[c].average_cov << '\t' << per_loc[c].min_cov << '\t' <<
			per_loc[c].max_cov << '\t' << per_loc[c].num_above_5x << '\t' << per_loc[c].num_above_10x << '\t' << per_loc[c].stack_depth << '\n';
	}
	avg_per_locus = avg_per_locus / num_loci;
	loc_sum.close();

	full_name = out_name + "_summary.txt";
	all_sum.open(full_name.c_str());
	all_sum << "Num individuals" << '\t' << per_ind.size() << '\n';
	all_sum << "Num loci" << '\t' << num_loci << '\n';
	all_sum << "Average per-locus coverage\t" << avg_per_locus << '\n';
	all_sum << "average per-individual coverage\t" << avg_per_ind << '\n';
	all_sum << "Num Loci with >=5x coverage \t" << num_above5x << '\n';
	all_sum << "Num loci with >=10x coverage \t" << num_above10x << '\n';
	all_sum << "Overall minimum per-locus coverage \t" << overall_min << '\n';
	all_sum << "Overall maximum per-locus coverage \t" << overall_max << '\n';
	all_sum.close();

	if (interactivemode)
	{
		cout << "\nDone! Enter Integer to quit: \n";
		cin >> end;
	}
	return 0;
}
