//
//  main.cpp
//  Lab2_Alignment
//
//  Created by 정구열 on 2018. 10. 17..
//  Copyright © 2018년 guyeol_jeong. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
class SequenceAlignment {
private:
  int match, mismatch, gap;
  int **gScore, **lScore;
  int **gBacktrack, **lBacktrack;
  int maxCol, maxRow, maxLScore;
  string v, w;
public:
  //Constructor for Sequence Alignment
  SequenceAlignment(string _seq1, string _seq2, int _match, int _mismatch, int _gap) {
    match = _match;
    mismatch = _mismatch;
    gap = _gap;
    v = _seq1;
    w = _seq2;
    LCSBackTrack();
  }
  ~SequenceAlignment() {
    for(int i = 0; i <= v.length(); ++i) {
      delete [] gScore[i];
      delete [] lScore[i];
      delete [] gBacktrack[i];
      delete [] lBacktrack[i];
    }
    delete [] gScore;
    delete [] lScore;
    delete [] gBacktrack;
    delete [] lBacktrack;
  }
  
  void LCSBackTrack() {
    gScore = new int*[v.length()+1];
    lScore = new int*[v.length()+1];
    gBacktrack = new int*[v.length()+1];
    lBacktrack = new int*[v.length()+1];
    for(int i = 0; i <= v.length(); ++i) {
      gScore[i] = new int[w.length()+1];
      lScore[i] = new int[w.length()+1];
      gBacktrack[i] = new int[w.length()+1];
      lBacktrack[i] = new int[w.length()+1];
    }
    gScore[0][0] = 0;
    lScore[0][0] = 0;
    maxLScore = lScore[0][0];
    maxCol = 0;
    maxRow = 0;
    for(int i = 1; i <= v.length(); i++) {gScore[i][0] = gap * i; lScore[i][0] = 0; gBacktrack[i][0] = 1;}
    for(int i = 1; i <= w.length(); i++) {gScore[0][i] = gap * i; lScore[0][i] = 0; gBacktrack[0][i] = 2;}
    for(int i = 1; i <= v.length(); i++) {
      for(int j = 1; j <= w.length(); j++) {
        if(v[i-1] == w[j-1]) {
          gScore[i][j] = max(gScore[i-1][j] + gap, max(gScore[i][j-1] + gap, gScore[i-1][j-1] + match));
          lScore[i][j] = max(lScore[i-1][j] + gap, max(lScore[i][j-1] + gap, max(lScore[i-1][j-1] + match, 0)));
        }
        else {
          gScore[i][j] = max(gScore[i-1][j] + gap, max(gScore[i][j-1] + gap, gScore[i-1][j-1] + mismatch));
          lScore[i][j] = max(lScore[i-1][j] + gap, max(lScore[i][j-1] + gap, max(lScore[i-1][j-1] + mismatch, 0)));
        }
        if(gScore[i][j] == gScore[i-1][j] + gap)
          gBacktrack[i][j] = 1;
        else if(gScore[i][j] == gScore[i][j-1] + gap)
          gBacktrack[i][j] = 2;
        else if(gScore[i][j] == gScore[i-1][j-1] + (v[i-1] == w[j-1]) ? match : mismatch)
          gBacktrack[i][j] = 3;
        if(lScore[i][j] == lScore[i-1][j] + gap){
          lBacktrack[i][j] = 1;
        } else if(lScore[i][j] == lScore[i][j-1] + gap) {
          lBacktrack[i][j] = 2;
        } else if(lScore[i][j] == lScore[i-1][j-1] + match) {
          lBacktrack[i][j] = 3;
        } else {
          lBacktrack[i][j] = 0;
        }
        if(lScore[i][j] > maxLScore) {maxLScore = lScore[i][j];maxCol = i;maxRow = j;}
      }
    }
  }
  
  void writeResult(string outPath){
    ofstream writeFile(outPath.data());
    int vl = (int)v.length();
    int wl = (int)w.length();
    int matchCnt = 0, misCnt = 0, gapCnt = 0;
    vector<char> vStack, wStack;
    while(vl && wl) {
      if(gBacktrack[vl][wl] == 1) {
        wStack.push_back('-');
        vStack.push_back(v[vl-1]);
        vl -= 1;
        gapCnt += 1;
      } else if(gBacktrack[vl][wl] == 2) {
        wStack.push_back(w[wl-1]);
        vStack.push_back('-');
        wl -= 1;
        gapCnt += 1;
      } else if(gBacktrack[vl][wl] == 3) {
        wStack.push_back(w[wl-1]);
        vStack.push_back(v[vl-1]);
        if(v[vl-1] == w[wl-1]) matchCnt += 1;
        else misCnt += 1;
        wl -= 1;
        vl -= 1;
      }
    }
    if(writeFile.is_open()) {
      writeFile << "1. Global Alignment\n";
      writeFile << "matches    : \t"<< matchCnt << "\n";
      writeFile << "mismatches : \t"<< misCnt << "\n";
      writeFile << "gaps       : \t"<< gapCnt << "\n";
      writeFile << "score      : \t"<< gScore[v.length()][w.length()] << "\n";
      writeFile <<"\n[Before Alignment]\n";
      writeFile << v << "\n";
      writeFile << w << "\n";
      writeFile << "\n[After Global Alignment]\n";
      // Print out 1st sequence
      while(!vStack.empty()) {
        writeFile << vStack.back();
        vStack.pop_back();
      }
      writeFile << "\n";
      while(!wStack.empty()) {
        writeFile << wStack.back();
        wStack.pop_back();
      }
      
      vl = maxCol;
      wl = maxRow;
      matchCnt = 0; misCnt = 0; gapCnt = 0;
      while(lBacktrack[vl][wl]!=0) {
        if(lBacktrack[vl][wl] == 1) {
          wStack.push_back('-');
          vStack.push_back(v[vl-1]);
          vl -= 1;
          gapCnt += 1;
        } else if(lBacktrack[vl][wl] == 2) {
          wStack.push_back(w[wl-1]);
          vStack.push_back('-');
          wl -= 1;
          gapCnt += 1;
        } else if(lBacktrack[vl][wl] == 3) {
          wStack.push_back(w[wl-1]);
          vStack.push_back(v[vl-1]);
          if(v[vl-1] == w[wl-1]) matchCnt += 1;
          else misCnt += 1;
          wl -= 1;
          vl -= 1;
        }
      }
      writeFile << "\n\n2. Local Alignment\n";
      writeFile << "matches    : \t"<< matchCnt << "\n";
      writeFile << "mismatches : \t"<< misCnt << "\n";
      writeFile << "gaps       : \t"<< gapCnt << "\n";
      writeFile << "score      : \t"<< maxLScore << "\n";
      writeFile <<"\n[Before Alignment]\n";
      writeFile << v << "\n";
      writeFile << w << "\n";
      writeFile << "\n[After Local Alignment]\n";
      while(!vStack.empty()) {
        writeFile << vStack.back();
        vStack.pop_back();
      }
      writeFile << "\n";
      while(!wStack.empty()) {
        writeFile << wStack.back();
        wStack.pop_back();
      }
      writeFile << "\n";
    }
    writeFile.close();
  }
};


int main(int argc, const char * argv[]) {
  string seq_file(argv[1]), seq_file2(argv[2]);
  string score_file(argv[3]), outputPath(argv[4]);
  string seq1, seq2;
  ifstream openFile(seq_file.data());
  ifstream openFile2(seq_file2.data());
  ifstream openFile3(score_file.data());
  
  if(openFile.is_open()) {
    string line;
    getline(openFile, line);
    getline(openFile, seq1);
    for(int i=0;i<=seq1.length();i++) {
      if(seq1[i]>=65 && seq1[i]<=92) seq1[i]+=32;
    }
    openFile.close();
  }
  if(openFile2.is_open()) {
    string line;
    getline(openFile2, line);
    getline(openFile2, seq2);
    for(int i=0;i<=seq2.length();i++) {
      if(seq2[i]>=65 && seq2[i]<=92) seq2[i]+=32;
    }
    openFile2.close();
  }
  if(openFile3.is_open()) {
    string line;
    getline(openFile3, line);
    int i = 0;
    while(line[i] != '=') i++;
    line.erase(0, i + 1);
    int match = stoi(line);
    getline(openFile3, line);
    i = 0;
    while(line[i] != '=') i++;
    line.erase(0, i + 1);
    int mismatch = stoi(line);
    getline(openFile3, line);
    i = 0;
    while(line[i] != '=') i++;
    line.erase(0, i + 1);
    int gap = stoi(line);
    openFile3.close();
    SequenceAlignment S(seq1, seq2, match, mismatch, gap);
    S.writeResult(outputPath);
  }
  return 0;
}
