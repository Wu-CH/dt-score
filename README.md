# dt-score

DT-Score (Delaunay Triangulation-Score)

## Introduction

DT-Score (Delaunay Triangulation-Score), aims to maximize the coverage of a given sensing area with obstacles. 
The DT-Score consists of two phases. In the first phase, we use a contour-based deployment to eliminate the coverage holes near the boundary of sensing area and obstacles. 
In the second phase, a deployment method based on the Delaunay Triangulation is applied for the uncovered regions.

## Build

```sh
git clone https://github.com/Wu-CH/dt-score.git
cd dt-score/source
make clean
make all
```

## Usage

- deployment

```sh
cd dt-score/test
../source/DT-Score/dt-score case_1_600
```

- calculage coverage

```sh
cd dt-score/test
../source/coverage/coverage case_1_600_result_DT-Score
```

## Citation

```txt
Chun-Hsien Wu, Kuo-Chuan Lee, Yeh-Ching Chung,
A Delaunay Triangulation based method for wireless sensor network deployment,
Computer Communications,
Volume 30, Issues 14–15,
2007,
Pages 2744-2752,
ISSN 0140-3664,
https://doi.org/10.1016/j.comcom.2007.05.017.
(http://www.sciencedirect.com/science/article/pii/S0140366407002046)
Keywords: Wireless sensor network; Sensor deployment; Sensor coverage; Obstacles; Delaunay Triangulation
```
