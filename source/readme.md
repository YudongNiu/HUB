# Hub Discovery from Hidden Networks in Knowledge Graphs

Compile
-------
Before compile the program, set the key parameters for our program in file "param.h".

**There are 5 parameters:**
* K: the size of the KMV sketch
* L: the number of sketch propagations
* LOOP: take the average result of LOOP executions (for approximate HUD queries)
* QNUM: the number of personalized queries sampled for each candidate meta-path (for personalized HUD queries)
* PATHCOUNT: the number of meta-paths sampled for scalability test with lengh 1, 2, 3 respectively

Then, compile the code for HUD problem by executing the following command on linux:

```sh
g++ -O3 main.py -o HUD
```

Running code
-------

To run the code for HUD problem, execute the following command on linux:

```sh
./HUD dataset centrality $\lambda$ $\beta$ len type 
```

**There are 6 parameters:**
* dataset: name of the KG
* centrality: take values of {d,df1,dp,h,hf1,hp}
* $\lambda$: parameter for HUD problem
* $\beta$: parameter for early termination
* len: run program on meta-paths of length
* type: the type of experiments to run


For example, the following command execute the sketch propagation framework to discover top-0.01 nodes in each candidate meta-path 

```sh
./HUD imdb df1 0.01 0 0 effect_prop_global_cross
```



Input Files
-----------
**The program HUD requires 4 input files:**
* node.dat stores nodes in KG. Each row in the file denotes an edge in the network.
* link.dat stores edges in KG. Row $i$ records the labels of node $v_i$.
* meta.dat stores the number of nodes in KG.
* dataset-cod-global-rules.dat stores the meta-paths mined by AnyBURL

**Additional input files:**
* dataset-cod-global-rules.dat1 stores meta-paths of length 1
* dataset-cod-global-rules.dat2 stores meta-paths of length 2
* dataset-cod-global-rules.dat3 stores meta-paths of length 3