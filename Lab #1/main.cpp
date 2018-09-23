#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <algorithm>

using namespace std;

/*
 The DNA sequences are each matched as follows.
 A = 0, C = 1, G = 2, T = 3
 */
class gibbsSampling {
private:
  int motifLen;
  int numOfString;
  int maxLength;
  int bestScore;
  int tmpScore;
  vector<string> dna;
  vector<string> bestMotifs;
  vector<vector<float> > profile;
public:
  gibbsSampling(int k, int t, int m, vector<string> _dna) {
    motifLen = k;
    numOfString = t;
    maxLength = m;
    dna = _dna;
    bestMotifs = randomlySelectMotif(dna, motifLen);
    bestScore = getScore(bestMotifs);
    for(int i = 0; i < 20; i++) {
      vector<string> motifs = gibbsSampler(dna);
      if(tmpScore < bestScore) {
        bestMotifs = motifs;
        bestScore = tmpScore;
      }
    }
    printOut();
  };
  
  vector<string> gibbsSampler(vector<string> _dna) {
    vector<string> motifs = randomlySelectMotif(_dna, motifLen);
    vector<string> _bestMotifs = motifs;
    int _bestScore = getScore(_bestMotifs);
    for(int i = 1; i <= numOfString; i++) {
      // Randomly extracts the index of the motif to exclude from the profile matrix.
      int j = (rand() % numOfString);
      profile = getProfileMatrix(motifs, j);
      string newMotif = _dna.at(j).substr(getMaxProfileMotif(_dna.at(j)), motifLen);
      motifs.at(j) = newMotif;
      int score = getScore(motifs);
      if(score < _bestScore) {
        _bestMotifs = motifs;
        _bestScore = score;
      }
    }
    tmpScore = _bestScore;
    return _bestMotifs;
  };
  
  vector<string> randomlySelectMotif(vector<string> _dna, int k) {
    vector<string> randomMotif;
    for(int i = 0; i < numOfString; i++) {
      int startIdx = rand() % (_dna.at(i).size() - k);
      string _randMotif = _dna.at(i).substr(startIdx, k);
      randomMotif.push_back(_randMotif);
    }
    return randomMotif;
  };
  
  int getScore(vector<string> motifs) {
    int score = 0;
    for(int i = 0; i < motifLen; i++) {
      int aSum = 0, cSum = 0, gSum = 0, tSum = 0;
      for(int j = 0; j < motifs.size(); j++) {
        if(motifs.at(j).at(i) == 'A' || motifs.at(j).at(i) == 'a') aSum ++;
        else if(motifs.at(j).at(i) == 'C' || motifs.at(j).at(i) == 'c') cSum ++;
        else if(motifs.at(j).at(i) == 'G' || motifs.at(j).at(i) == 'g') gSum ++;
        else if(motifs.at(j).at(i) == 'T' || motifs.at(j).at(i) == 't') tSum ++;
      }
      score += (aSum + cSum + gSum + tSum) - max(aSum, max(cSum, max(gSum, tSum)));
    }
    return score;
  };
  
  vector<vector <float> > getProfileMatrix(vector<string> motifs, int k) {
    vector<vector<float> > _profile(4, vector<float>(motifLen, 0.0));
    vector<vector <int> > count(4, vector<int>(motifLen, 1));
    for(int i = 0; i < motifLen; i++) {
      int total = 0;
      for(int j = 0; j < motifs.size(); j++) {
        if(j == k)
          continue;
        if(motifs.at(j).at(i) == 'A' || motifs.at(j).at(i) == 'a')
          total += ++count[0][i];
        else if(motifs.at(j).at(i) == 'C' || motifs.at(j).at(i) == 'c')
          total += ++count[1][i];
        else if(motifs.at(j).at(i) == 'G' || motifs.at(j).at(i) == 'g')
          total += ++count[2][i];
        else if(motifs.at(j).at(i) == 'T' || motifs.at(j).at(i) == 't')
          total += ++count[3][i];
      }
      for(int u = 0; u < 4; u++) _profile[u][i] = (float)count[u][i] / (float)(total + 4);
    }
    return _profile;
  };
  
  int getMaxProfileMotif(string dna) {
    float prof = 1.0;
    vector<float> _profileArr;
    string subStr;
    for(int i = 0; i < dna.length() - motifLen + 1; i++) {
      subStr = dna.substr(i, motifLen);
      for(int j = 0; j < motifLen; j++) {
        if(subStr.at(j) == 'A' || subStr.at(j) == 'a')
          prof *= profile[0][j];
        else if(subStr.at(j) == 'C' || subStr.at(j) == 'c')
          prof *= profile[1][j];
        else if(subStr.at(j) == 'G' || subStr.at(j) == 'g')
          prof *= profile[2][j];
        else if(subStr.at(j) == 'T' || subStr.at(j) == 't')
          prof *= profile[3][j];
      }
      _profileArr.push_back(prof);
    }
    float max = 0.0;
    int maxIdx = 0;
    for(int i = 0; i < _profileArr.size(); i++) {
      if(max < _profileArr.at(i)) {
        max = _profileArr.at(i);
        maxIdx = i;
      }
    }
    return maxIdx;
  }
  void printOut() {
    string filePath = "/Users/jeong-guyeol/Developer/BioInfomatics/Lab1_Motif_Finding/output.txt";
    ofstream writeFile(filePath.data());
    if( writeFile.is_open() ){
      for(vector<string>::iterator it = bestMotifs.begin(); it != bestMotifs.end(); it++)
        writeFile << *it << endl;
      for(vector<string>::iterator it = bestMotifs.begin(); it != bestMotifs.end(); it++) {
        for(int i = 0; i < motifLen; i++) {
          writeFile.setf(ios::fixed);
          writeFile.precision(6);
          if(it->at(i) == 'A' || it->at(i) == 'a')
            writeFile << profile[0][i] << '\t';
          else if(it->at(i) == 'C' || it->at(i) == 'c')
            writeFile << profile[1][i] << '\t';
          else if(it->at(i) == 'G' || it->at(i) == 'g')
            writeFile << profile[2][i] << '\t';
          else if(it->at(i) == 'T' || it->at(i) == 't')
            writeFile << profile[3][i] << '\t';
        }
        writeFile << endl;
      }
      writeFile.close();
    }
  };
};

int main(int argc, const char * argv[]) {
  string filePath = "/Users/jeong-guyeol/Developer/BioInfomatics/Lab1_Motif_Finding/test01.txt";
  ifstream openFile(filePath.data());
  vector<string> dna;
  if( openFile.is_open() ){
    string line;
    getline(openFile, line);
    int len = stoi(line);
    getline(openFile, line);
    int num = stoi(line);
    getline(openFile, line);
    int max = stoi(line);
    while(getline(openFile, line))
      dna.push_back(line);
    gibbsSampling gs(len, num, max, dna);
    openFile.close();
  }
  return 0;
}
