#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>

void getNames(std::vector <std::string> & indNames, const std::string & line);
void printVector(const std::vector <std::string> vecteur);
void treatLine(std::vector <std::string> & alignment, std::string line, const size_t nInd, const unsigned int cov);
void checkAlleles(std::string & refA, std::string & altA);
void addNucleotide(std::vector <std::string> & alignment, const std::string word, const size_t indPos, const std::string refA, const std::string altA, const size_t posDP, const size_t posAD, const size_t nFieldSep, const unsigned cov);
void writeFasta(const std::vector <std::string> & alignment, const std::vector <std::string> & indNames, const std::string outputName);

int main(int argc, char* argv[]){
	const std::string vcfFile(argv[1]); // name of the input vcf file
	const unsigned cov = atoi(argv[2]); // minimum number of reads to call an allele
	const std::string outputName(argv[3]); // minimum number of reads to call an allele
	
	std::size_t i(0);
	std::size_t nInd;
	std::string tmp_vector_string; // temporary vector of string to initialize 'alignment'
	std::vector <std::string> alignment; // contains the final output with sequences

	std::vector <std::string> indNames;
	std::ifstream infile(vcfFile.c_str());

	if(infile){ // if the vcf exists
		std::string line;

		while(std::getline(infile, line)){ // read the vcfFile
			if(line[0]=='#' & line[1] != '#'){ // read the header line with names of individuals
				getNames(indNames, line);
				nInd = indNames.size();

				for(i=0; i<2*nInd; ++i){
					alignment.push_back( tmp_vector_string );
				}
			}
			if(line[0]!='#'){ // after the header
				treatLine(alignment, line, nInd, cov);
			}
		} // end of: read the vcfFile
	}else{
		std::cout << "Error: the file " << vcfFile << " was not found" << std::endl;
		exit(EXIT_FAILURE);
	}

	writeFasta(alignment, indNames, outputName);

	return(0);
}

void getNames(std::vector <std::string> & indNames, const std::string & line){
	size_t i(0);
	std::string word;
	std::istringstream iss(line);
	while( std::getline( iss, word, '\t') ){
		++i;
		if(i > 9){
			indNames.push_back(word);
		}
	}
}

void printVector(const std::vector <std::string> vecteur){
	const size_t sizeV(vecteur.size());
	size_t i(0);

	for(i=0; i<sizeV; ++i){
		std::cout << "elements " << i << " : " << vecteur[i] <<  std::endl;
	}
}

void treatLine(std::vector < std::string> & alignment, std::string line, const size_t nInd, const unsigned int cov){
	std::size_t i(0); // 1-based counter
	std::size_t j(0); // 1-based counter
	std::size_t allele(0);

	std::size_t posDP(0); // 1-based position. If no DP is found in the FORMAT field --> remains to 0
	std::size_t posAD(0); // 1-based position. If no AD is found in the FORMAT field --> remains to 0
	
	std::size_t nFieldSep(0); // number of ':' in the field-9 (AD:DP:...)

	std::string contig("");
	std::string position("");
	std::string refA("");
	std::string altA("");

	std::string word;
	std::istringstream iss(line);
	while( std::getline( iss, word, '\t') ){
		++i;
		if(i == 1 ){
			contig = word;
		}
		
		if(i == 2 ){
			position = word;
		}

		if(i == 4 ){
			refA = word;
		}

		if(i == 5 ){
			altA = word;
			checkAlleles(refA, altA);
			if( altA[0] == '.'){ altA = refA ; }
		}

		if(i == 9 ){
			nFieldSep = std::count(word.begin(), word.end(), ':');

			std::istringstream iss2(word); // stream iss2 reads through word
			std::string word2;
			j = 0;
			
			while( std::getline( iss2, word2, ':')){
				++j;
				if(word2 == "DP"){ posDP=j; }
				if(word2 == "AD"){ posAD=j; }
			}
		}
		
		if(i >= 10){
			if( refA[0] == 'N' || altA[0] == 'N' ){ // if ref.size>1, or alt.size>1
				for(allele = 0; allele<2; ++allele){
					alignment[(i-10)*2 + allele].push_back('N');
				}
			}else{
				addNucleotide(alignment, word, i-10, refA, altA, posDP, posAD, nFieldSep, cov);
			}

		}
	}
}

void checkAlleles(std::string & refA, std::string & altA){
	// mask all vcf entries where variants are bigger than one nucleotide
	if( refA.size() == 1 ){
		if( altA.size() > 1 ){
			altA = 'N';
		}
	}else{
		refA = 'N';
	}
	
	if( refA[0] == 'N' ){ altA[0] = 'N'; }
}

void addNucleotide(std::vector <std::string> & alignment, const std::string word,  const size_t indPos, const std::string refA, const std::string altA, const size_t posDP, const size_t posAD, const size_t nFieldSep, const unsigned cov){
	// word = informations for one individual, for one position
	size_t i(0);
	size_t nFieldSep_tmp = std::count(word.begin(), word.end(), ':');

	unsigned int DP(0); 
	unsigned int cov_refA = 0;
	unsigned int cov_altA = 0;
	
	if( nFieldSep_tmp != nFieldSep ){
		alignment[indPos*2].push_back('N');
		alignment[indPos*2 + 1].push_back('N');
	}else{
		std::istringstream iss(word); // stream iss reads through word
		std::string word2;

		unsigned int test_altAll(0); // 0: no alternative allele; 1: possibly an alternative allele, upon coverage
				
		while( std::getline( iss, word2, ':')){
			++i;
			if( i==posDP ){
				if( word2[0] == '.' ){
					DP = 0;
				}else{
					DP = std::stoi(word2);
				}
			}
			
			if( i==posAD ){
				test_altAll = 1;
				size_t j(0);
				std::istringstream iss2(word2); // stream iss2 reads through word2
				std::string word3;
				while( std::getline( iss2, word3, ',')){
					++j;
					if( j==1 ){
						if( word3[0] == '.'){
							cov_refA = 0;
						}else{
							cov_refA = std::stoi(word3);
						}
					}
					
					if( j==2 ){
						if( word3[0] == '.'){
							cov_altA = 0;
						}else{
							cov_altA = std::stoi(word3);
						}
					}
				}
			}
		}
		if( test_altAll == 0 ){ // if only 'DP' in the FORMAT field (of VCF)
			if( DP < cov ){ // then homozygote 'N'
				alignment[indPos*2].push_back('N');
				alignment[indPos*2 + 1].push_back('N');
			}
			if( DP >= cov ){
				if( DP < 2*cov ){
					alignment[indPos*2].append(refA); // heterozygote 'ref'/'N'
					alignment[indPos*2 + 1].push_back('N');
				}else{
					alignment[indPos*2].append(refA); // homozygote 'ref'/'ref'
					alignment[indPos*2 + 1].append(refA);
				}
			}
		}
		if( test_altAll == 1){ // if 'AD' was also found in the FORMAT field (of VCF)
			if( cov_refA >= cov ){
				// if cov_refA >= cov
				alignment[indPos*2].append(refA); 
				if( cov_altA >= cov ){
					alignment[indPos*2 + 1].append(altA); 
				}else{
					if( cov_refA < 2*cov){
						alignment[indPos*2 + 1].push_back('N'); 
					}else{
						alignment[indPos*2 + 1].append(refA);
					}
				}
			}else{
				// if cov_refA < cov
				if( cov_altA < cov ){
					alignment[indPos*2].push_back('N');
					alignment[indPos*2 + 1].push_back('N');
				}else{
					alignment[indPos*2].append(altA);
					if( cov_altA < 2*cov ){
						alignment[indPos*2 + 1].push_back('N');
					}else{
						alignment[indPos*2 + 1].append(altA);
					}
				}
			}
		}
	}
}


void writeFasta(const std::vector <std::string> & alignment, const std::vector <std::string> & indNames, const std::string outputName){
	size_t i(0);
	const size_t nInd(indNames.size());
	std::ofstream fastaFile;
	fastaFile.open(outputName.c_str(), std::ios::out);
	
	for(i=0; i<nInd; ++i){
		// allele 1
		fastaFile << ">" << indNames[i] << "|allele1" << std::endl;
		fastaFile << alignment[i*2] << std::endl;
		
		// allele 2
		fastaFile << ">" << indNames[i] << "|allele2" << std::endl;
		fastaFile << alignment[i*2 + 1] << std::endl;
	}
	
	fastaFile.close();
}
