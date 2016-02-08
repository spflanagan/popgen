//Author: Sarah P. Flanagan
//Date: 31 March 2015
//Purpose: identify loci that are more differentiated between north-south on same coast than compared to similar latitude on another coast.
//Input: batch_x.pop1-pop2_fst.tsv from Stacks populations..well actually the fst and summary files from calculate_stacksCIs program

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

class fst_identity
{
public:
	int locus;
	string scaffold;
	int bp;

	fst_identity()
	{
		locus = int();
		scaffold = string();
		bp = int();
	}
};

class fst_record
{
public:
	fst_identity id;
	double fst;
	double smooth_fst;

	fst_record()
	{
		id = fst_identity();
		fst = double();
		smooth_fst = double();
	}
};

class summary_stats
{
public:
	double mean;
	double sd;
	double ci95;

	summary_stats()
	{
		mean = double();
		sd = double();
		ci95 = double();
	}
};

class grouped_fst_records
{
public:
	vector<fst_record> record;
	string name;

	grouped_fst_records()
	{
		record.push_back(fst_record());
		name = string();
	}
};

class grouped_group
{
public:
	vector<grouped_fst_records> group;
	vector<fst_identity> locus_id;

	grouped_group()
	{
		group.push_back(grouped_fst_records());
		locus_id.push_back(fst_identity());
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

void parse_fst_file(ifstream& file, string& filename, vector<fst_record>& fst)
{
	int count;
	string line;

	file.open(filename);
	FileTest(file, filename);
	count = 0;
	while (universal_getline(file, line))
	{
		if (!file.eof())
		{
			if (count != 0)
			{
				fst.push_back(fst_record());
				stringstream ss(line);
				ss >> fst[count - 1].id.locus >> fst[count - 1].id.scaffold >> fst[count - 1].id.bp >> fst[count - 1].fst >> fst[count - 1].smooth_fst;
			}
			count++;
		}
	}
	file.close();
	cout << filename << " read!\n";
}

void calculate_fst_diff(vector<fst_record>& fst1, vector<fst_record>& fst2, vector<fst_record>& diff)
{
	size_t t, tt;
	int count;

	count = 0;
	
	for (t = 0; t < fst1.size(); t++)
	{
		for (tt = 0; tt < fst2.size(); tt++)
		{
			//find matching loci and calculate differences
			if (fst1[t].id.bp == fst2[tt].id.bp && fst1[t].id.scaffold == fst2[tt].id.scaffold && fst1[t].id.locus == fst2[tt].id.locus)
			{
				//cout << "Locus " << fst1[t].id.locus << ", on " << fst1[t].id.scaffold << ", at basepair " << fst1[t].id.bp << " compared.\n";
				diff.push_back(fst_record());
				diff[count].id.bp = fst1[t].id.bp;
				diff[count].id.scaffold = fst1[t].id.scaffold;
				diff[count].id.locus = fst1[t].id.locus;
				diff[count].fst = fst1[t].fst - fst2[tt].fst;
				diff[count].smooth_fst = fst1[t].smooth_fst - fst2[tt].smooth_fst;
				count++;
			}
		}
	}
}

void calc_summary_stats(vector<fst_record>& diff, bool fst, summary_stats& diff_summ)
{
	int count;
	size_t t;

	if (fst == true)
	{
		count = 0;
		diff_summ.mean = 0;
		for (t = 0; t < diff.size(); t++)
		{
			diff_summ.mean = diff_summ.mean + diff[t].fst;
			count++;
		}
		diff_summ.mean = diff_summ.mean / count;
		diff_summ.sd = 0;
		for (t = 0; t < diff.size(); t++)
		{
			diff_summ.sd = diff_summ.sd + (diff[t].fst - diff_summ.mean)*(diff[t].fst - diff_summ.mean);
		}
		diff_summ.sd = sqrt(diff_summ.sd / count);
	}
	else
	{
		count = 0;
		diff_summ.mean = 0;
		for (t = 0; t < diff.size(); t++)
		{
			diff_summ.mean = diff_summ.mean + diff[t].smooth_fst;
			count++;
		}
		diff_summ.mean = diff_summ.mean / count;
		diff_summ.sd = 0;
		for (t = 0; t < diff.size(); t++)
		{
			diff_summ.sd = diff_summ.sd + (diff[t].smooth_fst - diff_summ.mean)*(diff[t].smooth_fst - diff_summ.mean);
		}
		diff_summ.sd = sqrt(diff_summ.sd / count);
	}
	diff_summ.ci95 = diff_summ.mean + 1.96*diff_summ.sd;
}

void compare_loci_to_reference(vector<fst_record>& diff, vector<fst_identity>& ref)
{
	size_t t, tt;
	bool found;

	for (t = 0; t < diff.size(); t++)
	{
		found = false;//assume this comparison isn't in the reference
		for (tt = 0; tt < ref.size(); tt++)
		{
			if (ref[tt].bp == diff[t].id.bp && ref[tt].scaffold == diff[t].id.scaffold && ref[tt].locus == diff[t].id.locus)
			{//if it matches, then it is in the reference
				found = true;
			}
		}
		if (found == false)
		{//if this comparison was not found in the reference, add it!
			ref.push_back(fst_identity());
			ref.back().bp = diff[t].id.bp;
			ref.back().scaffold = diff[t].id.scaffold;
			ref.back().locus = diff[t].id.locus;
		}
	}
}

int main()
{
	int i, end, count;
	size_t t, tt, ttt;
	vector<string> cline_file_names, north_file_names, south_file_names;
	string fst_difference_name, fst_diff_summary_name, sig_diff_name, cline_diff_name;
	vector<ifstream> cline_files, north_files, south_files;
	ofstream fst_difference, fst_diff_summary, sig_diff, cline_diff;
	vector<grouped_fst_records> cline_fsts, north_fsts, south_fsts;
	vector<grouped_fst_records> cline_to_cline, cline_to_north, cline_to_south;
	vector<summary_stats> summ_cc, summ_cn, summ_cs;
	int count_files_c, count_files_n, count_files_s, count_cc, count_cn, count_cs;
	string line;
	bool found;
	vector<bool> sig_locus;
	
	fst_difference_name = "fst_difference_values_4cline.txt";
	fst_diff_summary_name = "fst_difference_summary_4cline.txt";
	sig_diff_name = "fst_sig_diff_4cline.txt";
	cline_diff_name = "fsts_diff_in_4cline.txt";

	cline_file_names.push_back("E://ubuntushare//stacks//populations//fst_TXSP-ALST.txt");
	cline_file_names.push_back("E://ubuntushare//stacks//populations//fst_FLAB-FLCC.txt");
	cline_file_names.push_back("E://ubuntushare//stacks//populations//fst_FLSG-FLSI.txt");
	south_file_names.push_back("E://ubuntushare//stacks//populations//fst_TXSP-FLAB.txt");
	south_file_names.push_back("E://ubuntushare//stacks//populations//fst_TXSP-FLSI.txt");
	south_file_names.push_back("E://ubuntushare//stacks//populations//fst_FLAB-FLSI.txt");
	north_file_names.push_back("E://ubuntushare//stacks//populations//fst_FLCC-FLSG.txt");
	north_file_names.push_back("E://ubuntushare//stacks//populations//fst_ALST-FLSG.txt");
	north_file_names.push_back("E://ubuntushare//stacks//populations//fst_ALST-FLCC.txt");

	count_files_c = 0;
	for (t = 0; t < cline_file_names.size(); t++)
	{
		cline_files.push_back(ifstream());
		cline_fsts.push_back(grouped_fst_records());
		parse_fst_file(cline_files[t], cline_file_names[t], cline_fsts[t].record);
		cline_fsts[t].name = cline_file_names[t];
		count_files_c++;
	}
	count_files_n = 0;
	for (t = 0; t < north_file_names.size(); t++)
	{
		north_files.push_back(ifstream());
		north_fsts.push_back(grouped_fst_records());
		parse_fst_file(north_files[t], north_file_names[t], north_fsts[t].record);
		north_fsts[t].name = north_file_names[t];
		count_files_n++;
	}
	count_files_s = 0;
	for (t = 0; t < south_file_names.size(); t++)
	{
		south_files.push_back(ifstream());
		south_fsts.push_back(grouped_fst_records());
		parse_fst_file(south_files[t], south_file_names[t], south_fsts[t].record);
		south_fsts[t].name = south_file_names[t];
		count_files_s++;
	}
	
	count_cc = count_cn = count_cs = 0;
	for (t = 0; t < count_files_c; t++)
	{
		for (tt = t + 1; tt < count_files_c; tt++)
		{
			cline_to_cline.push_back(grouped_fst_records());
			summ_cc.push_back(summary_stats());
			cout << "Comparing " << cline_fsts[t].name << " and " << cline_fsts[tt].name << '\n';
			calculate_fst_diff(cline_fsts[t].record, cline_fsts[tt].record, cline_to_cline[count_cc].record);
			calc_summary_stats(cline_to_cline[count_cc].record, true, summ_cc[count_cc]);
			stringstream comp_name;
			comp_name << cline_fsts[t].name.substr(4, cline_fsts[t].name.size() - 4) << "_" << cline_fsts[tt].name.substr(4, cline_fsts[tt].name.size() - 4);
			cline_to_cline[count_cc].name = comp_name.str();
			count_cc++;
		}
		for (tt = 0; tt < count_files_n; tt++)
		{

			cout << "Comparing " << cline_fsts[t].name << " and " << north_fsts[tt].name << '\n';
			cline_to_north.push_back(grouped_fst_records());
			summ_cn.push_back(summary_stats());
			calculate_fst_diff(cline_fsts[t].record, north_fsts[tt].record, cline_to_north[count_cn].record);
			calc_summary_stats(cline_to_north[count_cn].record, true, summ_cn[count_cn]);
			stringstream comp_name;
			comp_name << cline_fsts[t].name.substr(4, cline_fsts[t].name.size() - 4) << "_" << north_fsts[tt].name.substr(4, north_fsts[tt].name.size() - 4);
			cline_to_north[count_cn].name = comp_name.str();
			count_cn++;
		}
		for (tt = 0; tt < count_files_n; tt++)
		{

			cout << "Comparing " << cline_fsts[t].name << " and " << south_fsts[tt].name << '\n';
			cline_to_south.push_back(grouped_fst_records());
			summ_cs.push_back(summary_stats());
			calculate_fst_diff(cline_fsts[t].record, south_fsts[tt].record, cline_to_south[count_cs].record);
			calc_summary_stats(cline_to_south[count_cs].record, true, summ_cs[count_cs]);
			stringstream comp_name;
			comp_name << cline_fsts[t].name.substr(4, cline_fsts[t].name.size() - 4) << "_" << south_fsts[tt].name.substr(4, south_fsts[tt].name.size() - 4);
			cline_to_south[count_cs].name = comp_name.str();
			count_cs++;
		}
	}
	cout << "\nConducted " << count_cc << " cline-cline comparisons, " << count_cn << " cline-north comparisons, and " << count_cs << " cline-south comparisons.\n";
	vector<fst_identity> all_loci;
	for (t = 0; t < cline_to_cline[0].record.size(); t++)
	{
		all_loci.push_back(fst_identity());
		all_loci[t].bp = cline_to_cline[0].record[t].id.bp;
		all_loci[t].scaffold = cline_to_cline[0].record[t].id.scaffold;
		all_loci[t].locus = cline_to_cline[0].record[t].id.locus;
	}
	for (t = 1; t < count_cc; t++)
		compare_loci_to_reference(cline_to_cline[t].record, all_loci);
	for (t = 1; t < count_cn; t++)
		compare_loci_to_reference(cline_to_north[t].record, all_loci);
	for (t = 1; t < count_cs; t++)
		compare_loci_to_reference(cline_to_south[t].record, all_loci);

	for (t = 0; t < all_loci.size(); t++)
		sig_locus.push_back(false);

	fst_difference.open(fst_difference_name);
	cout << "Writing to " << fst_difference_name << ".\n";
	fst_difference << "Locus\tScaffold\tBP";
	for (t = 0; t < count_cc; t++)
		fst_difference << '\t' << cline_to_cline[t].name.substr(0, cline_to_cline[t].name.size() - 4) << "_fst\t" << cline_to_cline[t].name.substr(0, cline_to_cline[t].name.size() - 4) << "_smooth";
	for (t = 0; t < count_cn; t++)
		fst_difference << '\t' << cline_to_north[t].name.substr(0, cline_to_north[t].name.size() - 4) << "_fst\t" << cline_to_north[t].name.substr(0, cline_to_north[t].name.size() - 4) << "_smooth";
	for (t = 0; t < count_cs; t++)
		fst_difference << '\t' << cline_to_south[t].name.substr(0, cline_to_south[t].name.size() - 4) << "_fst\t" << cline_to_south[t].name.substr(0, cline_to_south[t].name.size() - 4) << "_smooth";
	for (t = 0; t < all_loci.size(); t++)
	{
		fst_difference << '\n' << all_loci[t].locus << '\t' << all_loci[t].scaffold << '\t' << all_loci[t].bp;
		for (tt = 0; tt < count_cc; tt++)
		{
			found == false;
			for (ttt = 0; ttt < cline_to_cline[tt].record.size(); ttt++)
			{
				if (all_loci[t].bp == cline_to_cline[tt].record[ttt].id.bp && all_loci[t].locus == cline_to_cline[tt].record[ttt].id.locus && all_loci[t].scaffold == cline_to_cline[tt].record[ttt].id.scaffold)
				{
					fst_difference << '\t' << cline_to_cline[tt].record[ttt].fst << '\t' << cline_to_cline[tt].record[ttt].smooth_fst;
					//cout << "Locus " << all_loci[t].locus << ", on " << all_loci[t].scaffold << ", at basepair " << all_loci[t].bp << " compared.\n";
					found = true;
					if (cline_to_cline[tt].record[ttt].fst >= summ_cc[tt].ci95)
						sig_locus[t] = true;
				}
			}
			if (!found)
				fst_difference << '\t' << '\t';
		}

		for (tt = 0; tt < count_cn; tt++)
		{
			found == false;
			for (ttt = 0; ttt < cline_to_north[tt].record.size(); ttt++)
			{
				if (all_loci[t].bp == cline_to_north[tt].record[ttt].id.bp && all_loci[t].locus == cline_to_north[tt].record[ttt].id.locus && all_loci[t].scaffold == cline_to_north[tt].record[ttt].id.scaffold)
				{
					fst_difference << '\t' << cline_to_north[tt].record[ttt].fst << '\t' << cline_to_north[tt].record[ttt].smooth_fst;
					//cout << "Locus " << all_loci[t].locus << ", on " << all_loci[t].scaffold << ", at basepair " << all_loci[t].bp << " compared.\n";
					found = true;
					if (cline_to_north[tt].record[ttt].fst >= summ_cn[tt].ci95)
						sig_locus[t] = true;
				}
			}
			if (!found)
				fst_difference << '\t' << '\t';
		}

		for (tt = 0; tt < count_cs; tt++)
		{
			found == false;
			for (ttt = 0; ttt < cline_to_south[tt].record.size(); ttt++)
			{
				if (all_loci[t].bp == cline_to_south[tt].record[ttt].id.bp && all_loci[t].locus == cline_to_south[tt].record[ttt].id.locus && all_loci[t].scaffold == cline_to_south[tt].record[ttt].id.scaffold)
				{
					fst_difference << '\t' << cline_to_south[tt].record[ttt].fst << '\t' << cline_to_south[tt].record[ttt].smooth_fst;
					//cout << "Locus " << all_loci[t].locus << ", on " << all_loci[t].scaffold << ", at basepair " << all_loci[t].bp << " compared.\n";
					found = true;
					if (cline_to_south[tt].record[ttt].fst >= summ_cs[tt].ci95)
						sig_locus[t] = true;
				}
			}
			if (!found)
				fst_difference << '\t' << '\t';
		}
	}
	fst_difference.close();

	fst_diff_summary.open(fst_diff_summary_name);
	cout << "Writing to " << fst_diff_summary_name << ".\n";
	fst_diff_summary << "Comparison\tMean\tSD\t95CI";
	
	for (t = 0; t < count_cc; t++)
		fst_diff_summary << '\n' << cline_to_cline[t].name << '\t' << summ_cc[t].mean << '\t' << summ_cc[t].sd << '\t' << summ_cc[t].ci95;
	for (t = 0; t < count_cn; t++)
		fst_diff_summary << '\n' << cline_to_north[t].name << '\t' << summ_cn[t].mean << '\t' << summ_cn[t].sd << '\t' << summ_cn[t].ci95;
	for (t = 0; t < count_cs; t++)
		fst_diff_summary << '\n' << cline_to_south[t].name << '\t' << summ_cs[t].mean << '\t' << summ_cs[t].sd << '\t' << summ_cs[t].ci95;
	fst_diff_summary.close();
	
	sig_diff.open(sig_diff_name);
	cout << "Writing to " << sig_diff_name << ".\n";
	sig_diff << "Locus\tScaffold\tBP";
	for (t = 0; t < count_cc; t++)
		sig_diff << '\t' << cline_to_cline[t].name.substr(0, cline_to_cline[t].name.size() - 4) << "_fst\t" << cline_to_cline[t].name.substr(0, cline_to_cline[t].name.size() - 4) << "_smooth";
	for (t = 0; t < count_cn; t++)
		sig_diff << '\t' << cline_to_north[t].name.substr(0, cline_to_north[t].name.size() - 4) << "_fst\t" << cline_to_north[t].name.substr(0, cline_to_north[t].name.size() - 4) << "_smooth";
	for (t = 0; t < count_cs; t++)
		sig_diff << '\t' << cline_to_south[t].name.substr(0, cline_to_south[t].name.size() - 4) << "_fst\t" << cline_to_south[t].name.substr(0, cline_to_south[t].name.size() - 4) << "_smooth";
	for (t = 0; t < sig_locus.size(); t++)
	{
		if (sig_locus[t])
		{
			sig_diff << '\n' << all_loci[t].locus << '\t' << all_loci[t].scaffold << '\t' << all_loci[t].bp;

			for (tt = 0; tt < count_cc; tt++)
			{
				found = false;
				for (ttt = 0; ttt < cline_to_cline[tt].record.size(); ttt++)
				{
					if (all_loci[t].bp == cline_to_cline[tt].record[ttt].id.bp && all_loci[t].locus == cline_to_cline[tt].record[ttt].id.locus && all_loci[t].scaffold == cline_to_cline[tt].record[ttt].id.scaffold)
					{
						sig_diff << '\t' << cline_to_cline[tt].record[ttt].fst << '\t' << cline_to_cline[tt].record[ttt].smooth_fst;
						found = true;
					}
				}
				if (!found)
					sig_diff << '\t' << '\t';
			}
			found = false;
			for (tt = 0; tt < count_cn; tt++)
			{
				for (ttt = 0; ttt < cline_to_north[tt].record.size(); ttt++)
				{
					if (all_loci[t].bp == cline_to_north[tt].record[ttt].id.bp && all_loci[t].locus == cline_to_north[tt].record[ttt].id.locus && all_loci[t].scaffold == cline_to_north[tt].record[ttt].id.scaffold)
					{
						sig_diff << '\t' << cline_to_north[tt].record[ttt].fst << '\t' << cline_to_north[tt].record[ttt].smooth_fst;
						found = true;
					}
				}
				if (!found)
					sig_diff << '\t' << '\t';
			}
			found = false;
			for (tt = 0; tt < count_cs; tt++)
			{
				for (ttt = 0; ttt < cline_to_south[tt].record.size(); ttt++)
				{
					if (all_loci[t].bp == cline_to_south[tt].record[ttt].id.bp && all_loci[t].locus == cline_to_south[tt].record[ttt].id.locus && all_loci[t].scaffold == cline_to_south[tt].record[ttt].id.scaffold)
					{
						sig_diff << '\t' << cline_to_south[tt].record[ttt].fst << '\t' << cline_to_south[tt].record[ttt].smooth_fst;
						found = true;
					}
				}
				if (!found)
					sig_diff << '\t' << '\t';
			}
		}
	}
	sig_diff.close();

	//identify only those that are more different in cline-north or cline-south than cline-cline
	//grouped_group cline_different;
	//double this_cline_fst, this_comp_fst;

	//int num_comparisons = 0;
	//count = 0;
	//bool found_this_locus;
	//for (size_t v = 0; v < all_loci.size(); v++)
	//{
	//	for (t = 0; t < count_cc; t++)
	//	{
	//		cline_different.group.push_back(grouped_fst_records());
	//		cline_different.group[num_comparisons].name = cline_to_cline[t].name;
	//		num_comparisons++;
	//		found_this_locus = false;
	//		for (tt = 0; tt < cline_to_cline[t].record.size(); tt++)
	//		{
	//			if (all_loci[v].bp == cline_to_cline[t].record[tt].id.bp && all_loci[v].locus == cline_to_cline[t].record[tt].id.locus && all_loci[v].scaffold == cline_to_cline[t].record[tt].id.scaffold)
	//			{
	//				this_cline_fst = cline_to_cline[t].record[tt].fst;
	//				for (ttt = 0; ttt < count_cn; ttt++)
	//				{
	//					for (size_t iv = 0; iv < cline_to_north[ttt].record.size(); iv++)
	//					{
	//						if (cline_to_cline[t].record[tt].id.bp == cline_to_north[ttt].record[iv].id.bp && cline_to_cline[t].record[tt].id.locus == cline_to_north[ttt].record[iv].id.locus
	//							&& cline_to_cline[t].record[tt].id.scaffold == cline_to_north[ttt].record[iv].id.scaffold)
	//						{
	//							this_comp_fst = cline_to_north[ttt].record[iv].fst;
	//							if (abs(cline_to_north[ttt].record[iv].fst) > abs(cline_to_cline[t].record[tt].fst))
	//							{
	//								if (!found_this_locus)
	//								{
	//									found_this_locus = true;
	//									cline_different.locus_id.push_back(fst_identity());
	//									cline_different.locus_id[count].bp = cline_to_cline[t].record[tt].id.bp;
	//									cline_different.locus_id[count].scaffold = cline_to_cline[t].record[tt].id.scaffold;
	//									cline_different.locus_id[count].locus = cline_to_cline[t].record[tt].id.locus;
	//									cline_different.group[t].record[count].fst = cline_to_cline[t].record[tt].fst;
	//									cline_different.group[t].record[count].smooth_fst = cline_to_cline[t].record[tt].smooth_fst;
	//								}
	//								found = false;
	//								for (i = 0; i < num_comparisons; i++)
	//								{
	//									if (cline_different.group[i].name == cline_to_north[ttt].name)
	//										found = true;
	//								}
	//								if (!found)
	//								{
	//									cline_different.group.push_back(grouped_fst_records());
	//									num_comparisons++;
	//									cline_different.group[num_comparisons].name = cline_to_north[ttt].name;
	//								}
	//								cline_different.group[num_comparisons].record[count].fst = cline_to_north[ttt].record[iv].fst;
	//								cline_different.group[num_comparisons].record[count].smooth_fst = cline_to_north[ttt].record[iv].smooth_fst;
	//							}
	//							count++;
	//						}
	//					}
	//				}//for north
	//				for (ttt = 0; ttt < count_cs; ttt++)
	//				{
	//					for (size_t iv = 0; iv < cline_to_south[ttt].record.size(); iv++)
	//					{
	//						if (cline_to_cline[t].record[tt].id.bp == cline_to_south[ttt].record[iv].id.bp && cline_to_cline[t].record[tt].id.locus == cline_to_south[ttt].record[iv].id.locus
	//							&& cline_to_cline[t].record[tt].id.scaffold == cline_to_south[ttt].record[iv].id.scaffold)
	//						{
	//							this_comp_fst = cline_to_south[ttt].record[iv].fst;
	//							if (abs(cline_to_south[ttt].record[iv].fst) > abs(cline_to_cline[t].record[tt].fst))
	//							{
	//								if (!found_this_locus)
	//								{
	//									found_this_locus = true;
	//									cline_different.locus_id.push_back(fst_identity());
	//									cline_different.locus_id[count].bp = cline_to_cline[t].record[tt].id.bp;
	//									cline_different.locus_id[count].scaffold = cline_to_cline[t].record[tt].id.scaffold;
	//									cline_different.locus_id[count].locus = cline_to_cline[t].record[tt].id.locus;
	//									cline_different.group[t].record[count].fst = cline_to_cline[t].record[tt].fst;
	//									cline_different.group[t].record[count].smooth_fst = cline_to_cline[t].record[tt].smooth_fst;
	//								}
	//								found = false;
	//								for (i = 0; i < num_comparisons; i++)
	//								{
	//									if (cline_different.group[i].name == cline_to_south[ttt].name)
	//										found = true;
	//								}
	//								if (!found)
	//								{
	//									cline_different.group.push_back(grouped_fst_records());
	//									num_comparisons++;
	//									cline_different.group[num_comparisons].name = cline_to_south[ttt].name;
	//								}
	//								cline_different.group[num_comparisons].record[count].fst = cline_to_south[ttt].record[iv].fst;
	//								cline_different.group[num_comparisons].record[count].smooth_fst = cline_to_south[ttt].record[iv].smooth_fst;
	//							}
	//							count++;
	//						}
	//					}
	//				}//for south
	//			}
	//		}
	//	}
	//}
	//cline_diff.open(cline_diff_name);
	//cout << "Writing to " << cline_diff_name << '\n';

	//cline_diff.close();

	cout << "Done! Input integer to quit:\n";
	cin >> end;
	return 0;
}

