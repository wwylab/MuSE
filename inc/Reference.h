

#ifndef _REFERENCE_H_
#define _REFERENCE_H_

#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <boost/iostreams/device/mapped_file.hpp>
#include "htslib/sam.h"

class Reference {
	std::unordered_map<std::string, char*> chromosomes;
	std::unordered_map<std::string, long> chromosomesLen;
	std::unordered_map<std::string, std::string> chromosomesStr;
	boost::iostreams::mapped_file referenceMap;
	std::map<std::string, int> chrToIdx;
	std::vector<std::string> idxToChr;
	uint8_t* mNucleicAcidMap;
	char* mNucleicAcidMapStrictlyACGTN;
	bool ignoreNonACGT;
	bool keepOriginal;
	uint64_t refSize;

public:
	std::string fileName;

	Reference();
	~Reference();
	uint8_t convertToBamFormat(char in);
	uint8_t convertToStictlyACGTNFormat(char in);
	void setupChrSubstringBam(std::string *chrStr, const char *inputStr, int length);
	void setupChrSubstringStrictlyACGTN(std::string* chrStr, const char* inputStr, int length);
	void setupChrSubstring(std::string *chrStr, const char *inputStr, int length);
	void setupChrSubstringOriginal(std::string* chrStr, const char* inputStr, int length);
	void readFile(bool ifBamFormat);

	const std::unordered_map<std::string, std::string>& getAllChrStr() {return chromosomesStr;}
	const std::map<std::string, int>& getAllChrToIdx(){return chrToIdx;}
	uint64_t getRefSize(){return refSize;}
	uint32_t getIndexSize(){return chromosomesLen.size();}
	void setIgnoreNonACGT(bool _ignoreNonACGT) { ignoreNonACGT = _ignoreNonACGT; }
	void setKeepOriginal(bool _keepOriginal) {keepOriginal = _keepOriginal;}
	inline void openFile(std::string fileName_in) { referenceMap.open(fileName_in, boost::iostreams::mapped_file::readonly); fileName = fileName_in; }
	bool is_open(){return referenceMap.is_open();}
	const std::unordered_map<std::string, char*>& getAllChr(){return chromosomes;}
	
	inline char* getChr(const std::string& chrName) { 
		if (chromosomes.find(chrName) == chromosomes.end()){
			fprintf(stderr, "ERROR: Chr %s is not found in the reference\n", chrName.c_str());
       		exit(EXIT_FAILURE);
		}

		return chromosomes[chrName]; 
	}

	inline const char* getChrSubStr(const std::string& chrName, unsigned start){ 
		if (chromosomes.find(chrName) == chromosomes.end()){
			fprintf(stderr, "ERROR: Chr %s is not found in the reference\n", chrName.c_str());
       		exit(EXIT_FAILURE);
		}
		return chromosomes[chrName] + start;
	}

	inline long getChrLen(const std::string& chrName) { 
		if (chromosomesLen.find(chrName) == chromosomesLen.end()){
			fprintf(stderr, "ERROR: Chr %s is not found in the reference\n", chrName.c_str());
	   		exit(EXIT_FAILURE);
		}
		return chromosomesLen[chrName]; 
	}
	
	const std::string& getChrStr(const std::string& chrName) {
		if (chromosomesStr.find(chrName) == chromosomesStr.end()){
			fprintf(stderr, "ERROR: Chr %s is not found in the reference\n", chrName.c_str());
       		exit(EXIT_FAILURE);
		}
		return chromosomesStr[chrName];
	}
	
	inline char getChar(const std::string& chrName, unsigned start) {
		if (chromosomesStr.find(chrName) == chromosomesStr.end()){
			fprintf(stderr, "ERROR: Chr %s is not found in the reference\n", chrName.c_str());
       		exit(EXIT_FAILURE);
		}
		if (start >= chromosomesStr[chrName].size()){
			fprintf(stderr, "ERROR: Chr position %d is out of bound\n", start);
       		exit(EXIT_FAILURE);
		}

		return chromosomesStr[chrName][start];
	}

	int getIdxFromChr(const std::string& chrName) {
		auto iter = chrToIdx.find(chrName);
		if (iter == chrToIdx.end()){
			fprintf(stderr, "ERROR: Chr %s is not found in the reference\n", chrName.c_str());
       		exit(EXIT_FAILURE);
		}
		return iter->second;
	}

	const std::string& getChrName(int i){
		if (i >= idxToChr.size()){
			fprintf(stderr, "ERROR: Chr idx %d out of bound\n", i);
       		exit(EXIT_FAILURE);
		}
		return idxToChr[i];
	}

};

#endif