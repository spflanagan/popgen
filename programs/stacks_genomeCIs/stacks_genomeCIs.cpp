//Author: Sarah P. Flanagan, sflanagan@bio.tamu.edu
//Date: 22 August 2014
//Purpose: take output from stacks' populations program and generate genome-wide confidence intervals for the
//smoothed Fsts. This output can then be run through stacks_fst_plots.R and Fst graphs for each chromosome can be produced.
//Input file should be a batch_x.fst_Y-Z.tsv file. This is one of the output files from the populations module STACKS (Catchen et al. 2011)
//This code was developed using output from STACKS v.1.20, but should work with output from any of versions of the STACKS populations module
//***************************************
// Usage:
// -i:	Stacks Fst Input Filename (include path)
// -o:	path where output files will be written (include final backslash)
// -h: Help
// no arguments: Interactive
//***************************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>

using namespace std;

bool FileTest(ifstream& file, string filename, bool mode)
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
			mode = true;
			file.open(filename);
		}
	}
	return mode;

}

class fst_dat
{
public:
	int locus, bp;
	string chrom;
	double fst, smoothfst;

	//Class constructor
	fst_dat()
	{
		locus = 0;
		bp = 0;
		chrom = "";
		fst = 0;
		smoothfst = 0;
	}
};

class summ_dat
{
public:
	string chrom_name;
	double mean_fst, mean_smooth_fst, CI99_smooth, CI95_smooth, CI99_fst, CI95_fst, var_smooth, var_fst;
	int num_loci;

	//Class constructor
	summ_dat()
	{
		chrom_name = "";
		mean_fst = 0;
		mean_smooth_fst = 0;
		CI99_smooth = 0;
		CI95_smooth = 0;
		CI99_fst = 0;
		CI95_fst = 0;
		num_loci = 0;
		var_smooth = 0;
		var_fst = 0;
	}
};

int main(int argc, char* argv[])
{
	int end, count, line_count, z;
	string fst_input_name, fst_output_name, line, summ_output_name, out_path;
	ifstream fst_input;
	ofstream fst_output, summ_output;
	istringstream lin;
	int a, b, f, g;
	double h, i, j, k, l, m, n, o, p, q, r, s, t, v;
	string c, d, e, header;
	vector <fst_dat> data;
	summ_dat overall;
	bool interactivemode = false;
	

	string query;
	string tempstring1, tempstring2;

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nCalculate Genome-wide Confidence Intervals from Stacks Fst\n";
			cout << "-i:\tinput file (with path)\n";
			cout << "-o:\toutput path (with final backslash)\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nCalculate Genome-wide Confidence Intervals from Stacks Fst\n";
			cout << "-i:\tinput file (with path)\n";
			cout << "-o:\toutput path (with final backslash)\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	fst_input_name = "default";
	fst_output_name = "default";

	for (z = 1; z < argc - 1; z++)
	{
		tempstring1 = argv[z];
		tempstring2 = argv[z + 1];
		if (tempstring1 == "-i")
			fst_input_name = tempstring2;
		if (tempstring1 == "-o")
			out_path = tempstring2;
	}

	if (fst_input_name == "default")
	{
		cout << "Input Filename:\n";
		cin >> fst_input_name;
		interactivemode = true;
	}

	if (out_path == "default")
	{
		cout << "\nOutput Path:\n";
		cin >> out_path;
		interactivemode = true;
	}

	if (out_path.at(out_path.length() - 1) != '\\')
		out_path = out_path + "\\";

	
	cout << "\n\nInput File:\t" << fst_input_name;
	cout << "\nOutput Path:\t" << out_path;

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

	fst_input.open(fst_input_name);
	interactivemode = FileTest(fst_input, fst_input_name, interactivemode);
	line_count = 0;
	count = 0;
	while (getline(fst_input, line, '\n'))
	{
		lin.clear();
		lin.str(line);
		if (line_count == 0)
		{
			line_count++;
		}
		else
		{
			if (lin >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l >> m >> n >> o >> p >> q >> r >> s >> t >> v)//21
			{
				//batchID, locusID, Pop1id, pop2id, chrom, bp, col, overallpi, fst, fetp, oddsrat, cihigh, cilow, lod, 
				//corrfst, smoothfst, amovafst, corramovafst, smoothamovafst, smoothamovafstp, windowsnpcnt
				if (count == 0)
				{
					cout << "Reading fst file comparing pop " << c << " and pop " << d << ".\n";
					fst_output_name = out_path + "fst_" + c + "-" + d + ".txt";
					summ_output_name = out_path + "fst_" + c + "-" + d + "_summary.txt";
				}

				overall.mean_fst = overall.mean_fst + i;
				overall.mean_smooth_fst = overall.mean_smooth_fst + p;

				data.push_back(fst_dat());
				data[count].chrom = e;
				data[count].bp = f;
				data[count].locus = b;
				data[count].fst = i;
				data[count].smoothfst = p;
				count++;
				line_count++;
				
			}
			else
			{
				cout << "WARNING: failed to decode line " << line << ", count: " << count <<'\n';
			}
		}//end if count 0
	}//end while(getline)
	fst_input.close();

	cout << "File read. Calculating summary statistics and writing data to files.\n";
	overall.mean_fst = overall.mean_fst / count;
	overall.mean_smooth_fst = overall.mean_smooth_fst / count;
	overall.num_loci = count;

	fst_output.open(fst_output_name.c_str());
	fst_output << "Locus\tChrom\tBP\tFst\tSmoothFst\n";
	for (int i = 0; i < count; i++)
	{
		overall.var_fst = overall.var_fst + (overall.mean_fst - data[i].fst)*(overall.mean_fst - data[i].fst);
		overall.var_smooth = overall.var_smooth + (overall.mean_smooth_fst - data[i].smoothfst)*(overall.mean_smooth_fst - data[i].smoothfst);
		//write fst data to a file for R to read
		fst_output << data[i].locus << '\t' << data[i].chrom << '\t' << data[i].bp << '\t' << data[i].fst << '\t' << data[i].smoothfst << '\n';
	}
	fst_output.close();
	overall.var_fst = overall.var_fst / overall.num_loci;
	overall.var_smooth = overall.var_smooth / overall.num_loci;
	double std_dev_fst = sqrt(overall.var_fst);
	double std_dev_smooth = sqrt(overall.var_smooth);

	overall.CI99_fst = overall.mean_fst + 2.57583 * std_dev_fst;
	overall.CI95_fst = overall.mean_fst + 1.95996 * sqrt(overall.var_fst);
	overall.CI99_smooth = overall.mean_smooth_fst + 2.57583 * std_dev_smooth;
	overall.CI95_smooth = overall.mean_smooth_fst + 1.95996 * sqrt(overall.var_smooth);

	//write summary stats to a file
	summ_output.open(summ_output_name.c_str());
	summ_output << "LociNo\tMeanFst\tMeanSmoothFst\tCI95\tCI99\tCI95smooth\tCI99smooth\n";
	summ_output << overall.num_loci << '\t' << overall.mean_fst << '\t' << overall.mean_smooth_fst << '\t' << overall.CI95_fst << '\t'
		<< overall.CI99_fst << '\t' << overall.CI95_smooth << '\t' << overall.CI99_smooth << '\n';
	summ_output.close();

	if (interactivemode)
	{
		cout << "Done! Input integer to quit. ";
		cin >> end;
	}
	
	return 0;
}