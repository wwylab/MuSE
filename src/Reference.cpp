

#include "Reference.h"
#include <iostream>

using namespace std;

Reference::Reference():refSize(0) {
	mNucleicAcidMap = new uint8_t[122]();
	mNucleicAcidMap[int('A')] = 1;
	mNucleicAcidMap[int('C')] = 2;
	mNucleicAcidMap[int('G')] = 4;
	mNucleicAcidMap[int('T')] = 8;
	mNucleicAcidMap[int('N')] = 15;
	mNucleicAcidMap[int('a')] = 1;
	mNucleicAcidMap[int('c')] = 2;
	mNucleicAcidMap[int('g')] = 4;
	mNucleicAcidMap[int('t')] = 8;
	mNucleicAcidMap[int('n')] = 15;


	mNucleicAcidMap[int('R')] = 5;
	mNucleicAcidMap[int('Y')] = 10;
	mNucleicAcidMap[int('K')] = 12;
	mNucleicAcidMap[int('M')] = 3;
	mNucleicAcidMap[int('S')] = 6;
	mNucleicAcidMap[int('r')] = 5;
	mNucleicAcidMap[int('y')] = 10;
	mNucleicAcidMap[int('k')] = 12;
	mNucleicAcidMap[int('m')] = 3;
	mNucleicAcidMap[int('s')] = 6;

	mNucleicAcidMap[int('W')] = 9;
	mNucleicAcidMap[int('B')] = 14;
	mNucleicAcidMap[int('D')] = 13;
	mNucleicAcidMap[int('H')] = 11;
	mNucleicAcidMap[int('V')] = 7;
	mNucleicAcidMap[int('w')] = 9;
	mNucleicAcidMap[int('b')] = 14;
	mNucleicAcidMap[int('d')] = 13;
	mNucleicAcidMap[int('h')] = 11;
	mNucleicAcidMap[int('v')] = 7;

	ignoreNonACGT = false;
	keepOriginal = false;
	mNucleicAcidMapStrictlyACGTN = new char[122]();
	mNucleicAcidMapStrictlyACGTN[int('A')] = 'A';
	mNucleicAcidMapStrictlyACGTN[int('C')] = 'C';
	mNucleicAcidMapStrictlyACGTN[int('G')] = 'G';
	mNucleicAcidMapStrictlyACGTN[int('T')] = 'T';
	mNucleicAcidMapStrictlyACGTN[int('N')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('a')] = 'A';
	mNucleicAcidMapStrictlyACGTN[int('c')] = 'C';
	mNucleicAcidMapStrictlyACGTN[int('g')] = 'G';
	mNucleicAcidMapStrictlyACGTN[int('t')] = 'T';
	mNucleicAcidMapStrictlyACGTN[int('n')] = 'N';

	mNucleicAcidMapStrictlyACGTN[int('R')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('Y')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('K')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('M')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('S')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('r')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('y')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('k')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('m')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('s')] = 'N';

	mNucleicAcidMapStrictlyACGTN[int('W')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('B')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('D')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('H')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('V')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('w')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('b')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('d')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('h')] = 'N';
	mNucleicAcidMapStrictlyACGTN[int('v')] = 'N';
}

Reference::~Reference() {
	delete[] mNucleicAcidMap;
	delete[] mNucleicAcidMapStrictlyACGTN;
	for(auto& chr : chromosomes) {
		free(chr.second);
	}
	referenceMap.close();
}
uint8_t Reference::convertToBamFormat(char in) {
	return mNucleicAcidMap[in];
}

uint8_t Reference::convertToStictlyACGTNFormat(char in) {
	return mNucleicAcidMapStrictlyACGTN[in];
}

void Reference::setupChrSubstringBam(string* chrStr, const char* inputStr, int length) {
	for(int i = 0; i < length; i++) {
		chrStr->push_back(char(convertToBamFormat(inputStr[i])));
	}
}

void Reference::setupChrSubstringStrictlyACGTN(string* chrStr, const char* inputStr, int length) {
	for(int i = 0; i < length; i++) {
		chrStr->push_back(char(convertToStictlyACGTNFormat(inputStr[i])));
	}
}

void Reference::setupChrSubstring(string* chrStr, const char* inputStr, int length) {
	for(int i = 0; i < length; i++) {
		chrStr->push_back(char(toupper(inputStr[i])));
	}
}

void Reference::setupChrSubstringOriginal(string* chrStr, const char* inputStr, int length) {
	for(int i = 0; i < length; i++) {
		chrStr->push_back(char(inputStr[i]));
	}
}

void Reference::readFile(bool ifBamFormat) {
	const char* currPtr = referenceMap.const_data();
	const char* startChr = currPtr + 1;
	const char* endChr = static_cast<const char*>(rawmemchr(startChr, '\n'));
    const char* tempEndChr = static_cast<const char*>(memchr(startChr, '\t', endChr - startChr));
    if(tempEndChr == NULL) {
      tempEndChr = static_cast<const char*>(memchr(startChr, ' ', endChr - startChr));
    }
	if(tempEndChr != NULL) {
      endChr = tempEndChr;
	}
	string prevChr(startChr, endChr - startChr);
	
	int chrIdx = 0;
	currPtr = endChr;
	startChr = currPtr;
	endChr = static_cast<const char*>(rawmemchr(startChr, '\n'));
	currPtr = endChr + 1;

	while (currPtr - referenceMap.const_data() <= referenceMap.size()) {

		startChr = currPtr;
		endChr = static_cast<const char*>(memchr(startChr, '>', referenceMap.size() - (currPtr - referenceMap.const_data()) ));
		if (!endChr)
			endChr = referenceMap.const_data() + referenceMap.size();
		string* chrString = new string("");

		const char* newLinePtr = startChr;
		const char* lastLinePtr = startChr;

		do {
			newLinePtr = static_cast<const char*>(memchr(newLinePtr, '\n', endChr - newLinePtr));
			if (!newLinePtr)
				newLinePtr = endChr;

			if (ifBamFormat)
				setupChrSubstringBam(chrString, lastLinePtr, newLinePtr - lastLinePtr);
			else if (ignoreNonACGT)
				setupChrSubstringStrictlyACGTN(chrString, lastLinePtr, newLinePtr - lastLinePtr);
			else if (keepOriginal)
				setupChrSubstringOriginal(chrString, lastLinePtr, newLinePtr - lastLinePtr);
			else
				setupChrSubstring(chrString, lastLinePtr, newLinePtr - lastLinePtr);

			lastLinePtr = ++newLinePtr;
		} while(newLinePtr < endChr);
		
		chromosomesLen[prevChr] = chrString->size();
		refSize += chrString->size();
		chrToIdx.insert(make_pair(prevChr, chrIdx++));
		idxToChr.push_back(prevChr);
		if (ifBamFormat){
			char* chrStr = strdup(chrString->c_str());
			chromosomes[prevChr] = chrStr;
		}else
			chromosomesStr.insert(make_pair(prevChr, move(*chrString)));

		delete chrString;
		if(endChr == referenceMap.const_data() + referenceMap.size()) {
			return;
		}
		currPtr = endChr + 1;

		startChr = currPtr;
		endChr = static_cast<const char*>(rawmemchr(startChr, '\n'));
		tempEndChr = static_cast<const char*>(memchr(startChr, '\t', endChr - startChr));
		if(tempEndChr == NULL) {
			tempEndChr = static_cast<const char*>(memchr(startChr, ' ', endChr - startChr));
		}
		if(tempEndChr != NULL) {
			endChr = tempEndChr;
		}

		string chr(startChr, endChr - startChr);
		prevChr = chr;

		currPtr = endChr;
		startChr = currPtr;
		endChr = static_cast<const char*>(rawmemchr(startChr, '\n'));
		currPtr = endChr + 1;
	}

}
