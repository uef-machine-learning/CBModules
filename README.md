# CBModules
Implementation of the following algorithms:
* K-Means initialization (9 different methods) [1]
* Random Swap algorithm [2]
* Genetic Algorithm [3]

This package is also needed as a dependency to be able to build some of the projects in the UEF organization (like [Random Swap](https://github.com/uef-machine-learning/RandomSwap))

# Compile and Run

## K-means

**Compile:**
`make cbkm`

**Example with 15 clusters:**
```
./cbkm -O -S15 datasets/s1.txt output.txt
./cbkm -O -S15 datasets/s1.ts output.ts
```

**Steinley initialization method, two repeats:**
```
./cbkm -I5 -O -R2 -S15 datasets/s1.txt output.txt
```

**K-means++ initialization method, two repeats:**
`./cbkm -I2 -O -R2 -S15 datasets/s1.txt output.txt`


**Run ./cbkm for help:**
```
CBKM	Version 0.65	4.4.2017
Repeated K-means algorithm.
Usage: CBKM [-option] <training set> [initial cb/pa] <result codebook>
For example: CBKM bridge initial tmp

  Options:
  -Bn              = Sample percentage 10000/(1..10000) (0..10000, default=1000)
  -Hn1,n2          = 
                      n1: Hybrid methods combining K-means
                            0 = None (default)
                            1 = Kmeans+DensityPeaks
                      n2: Hybrid method variant (0..10000, default=1)
  -I[n1,n2,..n4]   = 
                      n1: Initialization method
                            0 = DensityPeaks
                            1 = Random (default)
                            2 = K-means++
                            3 = Bradley
                            4 = Projection
                            5 = Luxburg
                            6 = MaxMin
                            7 = Random partition (Forgy)
                            8 = Splitting
                            9 = Sorting heuristic
                      n2: Initialization method variant (Maxmin-variant: 1=rand point first, 2=mean first) (0..10000, default=1)
                      n3: Initialization method option (0..10000, default=0)
                      n4: Initialization method (third) option (0..10000, default=0)
  -O               = Overwrite existing file (default=NO)
  -P               = Save partition file (default=NO)
  -Qn              = Quiet level (0..99, default=2)
  -Rn              = Number of repeats (0..1000000, default=1)
  -Sn              = Number of clusters (1..5000, default=256)
  -Tn              = Iterations, 0=INF (0..50000000, default=0)
  -Zn              = Random seed (0=clock)  (0..2147483647, default=0)
```

## Random Swap

**Compile:**
`make cbrs`

**Example with 15 clusters:**
`./cbrs -O -S15 datasets/s1.txt output.txt`

Run `./cbrs`  for help

## Genetic Algorithm

**Compile:**
`make cbga`

**Example with 15 clusters:**
`./cbga -O -S15 datasets/s1.txt output.txt`

Run `./cbga` for help

**Debug**
Compile with debug flags:
`make cbkm DEBUG=-g`


# References
[1] P Fränti, S Sieranoja, "How much can k-means be improved by using better initialization
and repeats?", Pattern Recognition, 93, 95-112, 2019. https://doi.org/10.1016/j.patcog.2019.04.014

[2] P. Fränti, "Efficiency of random swap clustering", Journal of Big Data, 5:13, 1-29, 2018. (pdf) JF=1 https://doi.org/10.1186/s40537-018-0122-y

[3] P. Fränti, "Genetic algorithm with deterministic crossover for vector quantization", Pattern Recognition Letters, 21 (1), 61-68, 2000. (pdf) IF=1.03 https://doi.org/10.1016/S0167-8655(99)00133-6
