지능형 생물정보학 Lab #2
===================
### 이름 : Guyeol Jeong  학번 : 2012004236

## 1. 프로그램 설명
### 1.1 입력
* 첫 번째 인자로 첫 번째 Sequence 파일의 경로를 입력한다.
* 두 번째 인자로 두 번째 Sequence 파일의 경로를 입력한다.
* 세 번째 인자로 score 파일의 경로를 입력한다.
* 네 번째 인자로 결과가 쓰여질 output 파일의 경로를 입력한다.
1. Sequence 파일은
	* 한 줄로 된 최대 100글자 sequence 하나로 이루어져있다.
2. score 파일은
  * match, mismatch, gap에 따른 score를 순서대로 적어놓는다.
### 1.2 출력
1. Global Alignment에 따른 score와 sequence를 출력시킨다.
2. Local Alignment에 따른 score와 sequence를 출력시킨다.
## 2. 로직
### 두 가지 알고리즘 for global alignment (Needleman-Wunsch), for local alignment (Smith-Waterman)을 구현하였다.
* 알고리즘에 대한 자세한 설명은 밑에 링크를 보자
  * [Needleman-Wunsch ](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
  * [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
