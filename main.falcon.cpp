#include <iostream>
#include <fstream>
#include <cstring>
#include <string.h>
#include <vector>
#include <math.h>
#include <numeric>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <map>
#include <time.h>
#include <ctime>
#include <numeric>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <stdlib.h>
#include <iomanip>

using namespace std;
#define malloc2D(name, xDim, yDim, type) do {               \
    name = (type **)malloc(xDim * sizeof(type *));          \
    assert(name != NULL);                                   \
    name[0] = (type *)malloc(xDim * yDim * sizeof(type));   \
    assert(name[0] != NULL);                                \
    for (size_t i = 1; i < xDim; i++)                       \
        name[i] = name[i-1] + yDim;                         \
} while (0)

const int fval = 4;
const int max_vec_size = 596;//9;
using namespace std;

int getdir (string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        string filename=dirp->d_name;
        size_t found=filename.find(".dta");
        if(found!=string::npos){
        files.push_back(string(dirp->d_name));
        }
        }
    closedir(dp);
    return 0;
}

void h_gen_kmers(const vector<int> spectrum, int pos, int k, vector<int>& d_spectra)
{
	vector<int> kmer(spectrum.begin()+pos, spectrum.begin()+k+pos);
	for (int i=0; i<kmer.size(); ++i) 
	{
        //cout<< kmer[i]<< " ";
		d_spectra.push_back(kmer[i]);
    }
	//cout<<endl;

	return;
}



//void fset(std::map<string,vector<int>> &input, int k, int number_of_files);
void fset(int* input, int k, int total, int number_of_files, int* output);
void kmeans(int **spectra,      /* in: [numObjs][numCoords] */
			int *precursorMass, //in: [numObjs]
            int  epf,    /* entries per file numCords*/
            int  numSpectra,      /* no. of spectra  numObjs*/
            int  numClusters,  /* no. clusters */
            float threshold,    /* % spectra change membership */
            int   *membership,   /* out: [numSpectra] */
            int   loop_iterations);

//
//
//void computeHistogram(vector<int> input);
//
//
//cudaError_t addWithCuda(int *c, const int *a, const int *b, size_t size);


std::ostream& operator<<(std::ostream& o, const pair<float,string>& p)
{
  return o << "Filename: " << p.second << "\t" << "Precursor ion mass: " << p.first;
}


std::vector<int> createHistogram (std::vector<int> spectrum, float range, int bins)
{
	const float binSize = range/bins;
	std:vector<int> Histogram(bins);

	for (int i=0; i< spectrum.size(); i++)
	{
		int bin = int(floor(spectrum[i]/binSize));
		Histogram[bin] += 1;
	}

	return Histogram;
}


void fset(int* inputData, int* inputHistogram, int* spectraSizes, int totalHistogramElements, int totalElements, int totalSpectra, int fsetSize, 
			int blocks, int threads, int* membershipArray);

void printDeviceArray(int* input, int size);


static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    
	if (x == (n-1))
	{
		float ratio = 1;
		int   c      =  ratio * w;
 
		cout << setw(3) << (int)(ratio*100) << "% [";
		for (int x=0; x<c; x++) cout << "=";
		for (int x=c; x<w; x++) cout << " ";
		cout << "]\r" << flush;

		return;
	}
	
	else
	{
		//if ( (x != n) && (x % (n/100+1) != 0) )
			//return;
 
		float ratio  =  x/(float)n;
	
		int   c      =  ratio * w;
 
		cout << setw(3) << (int)(ratio*100) << "% [";
		for (int x=0; x<c; x++) cout << "=";
		for (int x=c; x<w; x++) cout << " ";
		cout << "]\r" << flush;

		return;
	}
}

int main(int argc, char* argv[]){

	if(argc!=3)
	{
		cout<<"Usage is ./aout [directoryname] [number of blocks]"<<endl;
		exit(0);
	}


	time_t start, end;
	double duration ;
	time(&start);

	cout << "This program is going to score the spectral data, specifically the fset scores pairwise\n";


	// Read filenames from disk
	string directory=string(argv[1]);
	vector<string> filenames = vector<string>();
	getdir(directory,filenames);

	// Read spectral data into vector containers
	typedef std::multimap<pair< float,string>, std::vector<int> > Spectrum; // pair(spectrum name, precursor ion mass), spectral data
	vector<int> spectraSizes;
	Spectrum rawData;
	
	// Reading files from disks
	std::cout << "Reading files from disk...\n";
	for (unsigned int i = 0;i <filenames.size();i++) 
	{
		float precursorMass;
		vector<int> mz;
		string currentFile = (directory + filenames[i]);
		//cout<< currentFile<<endl;
		ifstream file;
		string s;
		file.open(currentFile.c_str());
		int count = 0;
		while(!file.eof())
		{
			getline(file,s);
			

			if (count != 0)
			{
				//cout<<s<<endl;
				char * cstr = new char [s.length()+1];
				std::strcpy (cstr, s.c_str());
				char * p = std::strtok (cstr," ");

				if (p!=0)
				  {
					int pp = int (atof(p) +0.5);
					//cout << pp << '\n';
					mz.push_back(pp);
					//mz[count-1] = pp;
					p = std::strtok(NULL," ");
				  }

			}
			else
			{
				char * cstr = new char [s.length()+1];
				std::strcpy (cstr, s.c_str());
				char * p = std::strtok (cstr," ");

				precursorMass = (float(atof(p)));
			}
			count++;

			
			
		}
		//spectraSizes.push_back(mz.size());
		rawData.insert(make_pair(make_pair(precursorMass,filenames[i]),mz));
		
		loadbar(i,filenames.size());
	}

	//std::cout << "Spectral data is :\n";
	Spectrum::iterator it;
	int j = 0;
	for (it=rawData.begin(); it!=rawData.end(); ++it)
	{
	    //std::cout << (*it).first << " => " << endl;
		//for (int i=0; i<(*it).second.size(); i++)
			//cout << (*it).second[i] << '\t';
		//cout<<endl;
		//cout<< "Size of spectrum in peaks: "<< (*it).second.size() <<"\n";
		vector<int>::const_iterator it2;
		it2 = std::max_element((*it).second.begin(), (*it).second.end());
		  //cout << " the max is " << *it2 << endl;
		spectraSizes.push_back((*it).second.size());
	}
	cout <<endl;

	//Sort Spectra Data according to precursor ion masses
	//In current multimap format, the data is already sorted. Need to update this if change in implementation is done later on

	//Create Histogram for spectral data
	Spectrum histogram;
	const int bins = 500;
	const float range = 3000;

	for (it=rawData.begin(); it!=rawData.end(); ++it)
	{
		histogram.insert(make_pair( (*it).first , createHistogram((*it).second, range, bins) ));
	}

	//std::cout << "Corresponding histogram data is :\n";
	//for (it=histogram.begin(); it!=histogram.end(); ++it)
	{
		//std::cout << (*it).first << " => " << endl;
		//for (int i=0; i<(*it).second.size(); i++)
			//cout << (*it).second[i] << '\t';
		//cout<<endl;
	}


	// Transfer data to linear array for manipulation on the GPU
	int totalElements = std::accumulate(spectraSizes.begin(),spectraSizes.end(),0);
	int totalHistogramElements = bins*filenames.size();
	int *h_data = new int[totalElements];
	int *h_histogram = new int[totalHistogramElements];
	int sizeCounter = 0;
	int histogramCounter = 0;
	for (it=rawData.begin(); it!=rawData.end(); ++it)
	{
		for (int i=0; i<(*it).second.size(); i++)
			h_data[sizeCounter++] = (*it).second[i];
	}

	for (it=histogram.begin(); it!=histogram.end(); ++it)
	{
		for (int i=0; i<(*it).second.size(); i++)
			h_histogram[histogramCounter++] = (*it).second[i];
		
	}

	int* spectraSizesArray = new int[spectraSizes.size()];
	for (int i=0 ; i < spectraSizes.size(); i++)
		spectraSizesArray[i] = spectraSizes[i];

	int* scores = new int[filenames.size() * filenames.size()];

	cout << "Total number of spectra is: "<< filenames.size() <<endl;
	int totalSpectra = filenames.size();
	int blocks = atoi(argv[2]);
	int threads = totalSpectra/blocks;
	
	fset(h_data,h_histogram,spectraSizesArray, totalHistogramElements,totalElements, totalSpectra, fval,blocks, threads,  scores);
	


    return 0;
}
