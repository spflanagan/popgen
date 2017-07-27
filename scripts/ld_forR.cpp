#include <Rcpp.h>
#include <string>
#include <vector>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

class ld_info
{
public: 
  double d, dprime;
  string chrom1, chrom2, locus1, locus2;
  double pA, pB, pAB, paB, pAb, pab;
  
  ld_info()
  {
    d = dprime = double();
    chrom1 = chrom2 = locus1 = locus2 = string();
    pA = pB = pAB = paB = pAb = pab = double();
  }
};

class locus_info
{
public:
  string chrom, locusID;
  int bp;
  vector<int> maternal;
  vector<int> paternal;
  
  locus_info()
  {
    chrom = locusID = string();
    bp = int();
    maternal = paternal = vector<int>();
  }
};


ld_info AdultPopLD(locus_info & locus1, locus_info &locus2)
{
  //To do this, need to compare allele frequencies. 
  //D=x11-p1q1
  //x11 = observed freq of major allele in locus A and major allele in locus B
  //p1 = observed freq of major allele in locus A
  //q1 = observed freq of major allele in locus B
  //if D > 0, Dmax = min(p1q1, p2q2)
  //if D < 0, Dmax = min(p1q2, p2q1)
  //D' = D/Dmax
  ld_info result;
  int NumAlleles = 2;
  int num_adults, nad1, nad2;
  vector <int> num_allele_1;
  vector <int> num_allele_2;
  vector <double> freq_allele_1;
  vector <double> freq_allele_2;
  double joint_alleles[2][2];
  double D[2][2]; 
  double Dmax[2][2];
  int count_a = 0;
  int f, ff, fff, al, count;
  double Dprime, d_allele_avgs;
  int numA1B1 = 0;
  
  num_adults = nad1 = nad2 = 0;
 
  //figure out which loci are polymorphic
  for (al = 0; al < NumAlleles; al++)
  {
    num_allele_1.push_back(0);
    num_allele_2.push_back(0);
    freq_allele_1.push_back(0);
    freq_allele_2.push_back(0);
    for (f = 0; f < NumAlleles; f++)
    {
      joint_alleles[al][f] = 0;
      D[al][f] = 0;
      Dmax[al][f] = 0;
    }
  }
  
  int maternal_1, maternal_2, paternal_1, paternal_2;
  
  for (fff = 0; fff < locus1.maternal.size(); fff++)
  {
    maternal_1 = locus1.maternal[fff];
    maternal_2 = locus2.maternal[fff];
    paternal_1 = locus1.paternal[fff];
    paternal_2 = locus2.paternal[fff];
  
    if (maternal_1 >= 0 && maternal_2 >= 0)
    {			
      num_allele_1[maternal_1]++;
      num_allele_1[paternal_1]++;
      nad1++;
      num_allele_2[maternal_2]++;
      num_allele_2[paternal_2]++;
      nad2++; 
      joint_alleles[maternal_1][maternal_2]++;
      joint_alleles[paternal_1][paternal_2]++;
      num_adults++;
    }
    
  }//end for fff < PopulationSize
  
  int major_allele_1 = 0;
  int major_allele_2 = 0;
  for (al = 0; al < NumAlleles; al++)
  {
    freq_allele_1[al] = (num_allele_1[al]) / (2 * double(nad1));
    freq_allele_2[al] = (num_allele_2[al]) / (2 * double(nad2));
    if (freq_allele_1[al] > freq_allele_1[major_allele_1])
      major_allele_1 = al;
    if (freq_allele_2[al] > freq_allele_2[major_allele_2])
      major_allele_2 = al;
    for (f = 0; f < NumAlleles; f++)
      joint_alleles[al][f] = joint_alleles[al][f] / (2 * double(num_adults));
  }
  
  if (freq_allele_1[major_allele_1] != 1 && freq_allele_2[major_allele_1] != 1)//if it's polymorphic
  {
    d_allele_avgs = 0;
    count = 0;
    for (f = 0; f < NumAlleles; f++)
    {
      for (ff = 0; ff < NumAlleles; ff++)
      {
        if (freq_allele_1[f] > 0 && freq_allele_2[ff] > 0)
        {
          D[f][ff] = joint_alleles[f][ff] - freq_allele_1[f] * freq_allele_2[ff];
          d_allele_avgs = d_allele_avgs + fabs(D[f][ff]);
          count++;
          if (D[f][ff] < 0)
            Dmax[f][ff] = min(freq_allele_1[f] * freq_allele_2[ff], (1 - freq_allele_1[f])*(1 - freq_allele_2[ff]));
          else
            Dmax[f][ff] = min((1 - freq_allele_1[f])*freq_allele_2[ff], freq_allele_1[f] * (1 - freq_allele_2[ff]));
        }
      }
    }
    double dcount = count;
    d_allele_avgs = d_allele_avgs / dcount;
    
    Dprime = 0;
    bool decentDmax = true;
    for (f = 0; f < NumAlleles; f++)
    {
      for (ff = 0; ff < NumAlleles; ff++)
      {
        if (freq_allele_1[f] > 0 && freq_allele_2[ff] > 0)
        {
          if (Dmax[f][ff] > 0)
            Dprime = Dprime + freq_allele_1[f] * freq_allele_2[ff] * fabs(D[f][ff]) / Dmax[f][ff];
          else
            decentDmax = false;
          
        }
      }
    }
    if (!decentDmax)
      Dprime = -5;
    
  }
  
  if (Dprime > 1.1)
    cout << "D' greater than 1.";
  result.d = d_allele_avgs;
  result.dprime = Dprime;
  result.chrom1 = locus1.chrom;
  result.locus1 = locus1.locusID;
  result.chrom2 = locus2.chrom;
  result.locus2 = locus2.locusID;
  
  result.pA = freq_allele_1[0];
  result.pB = freq_allele_2[0];
  result.pAB = joint_alleles[0][0];
  result.pAb = joint_alleles[0][1];
  result.paB = joint_alleles[1][0];
  result.pab = joint_alleles[0][0];
  
  return result;
}//end Adult Pop LD

//MODULES
RCPP_MODULE(ld_module) {
  class_<ld_info>("ld_info")
  
  .constructor()
  
  .field("d", &ld_info::d)
  .field("dprime", &ld_info::dprime)
  .field("chrom1", &ld_info::chrom1)
  .field("chrom2", &ld_info::chrom2)
  .field("locus1", &ld_info::locus1)
  .field("locus2", &ld_info::locus2)
  .field("pA", &ld_info::pA)
  .field("pB", &ld_info::pB)
  .field("pAB", &ld_info::pAB)
  .field("paB", &ld_info::paB)
  .field("pAb", &ld_info::pAb)
  .field("pab", &ld_info::pab)
}

RCPP_MODULE(loc_module) {
  class_<locus_info>("locus_info")
  
  .constructor()
  
  .field("chrom", &loc_info::chrom)
  .field("locusID", &loc_info::locusID)
  .field("bp", &loc_info::bp)
  .field("maternal", &loc_info::maternal)
  .field("paternal", &loc_info::paternal)

}

RCPP_MODULE(calc_ld){
  using namespace Rcpp;
  function("AdultPopLD",&AdultPopLD,List::create(_["locus1"],_["locus2"]),
           "Calculates pairwise LD")
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
