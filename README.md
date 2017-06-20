# The Co-Evolution Model for Social Network Evolving and Opinion Migration

### Graph generation
Located in directory *./src/generation*.

### Inference 
Located in directory *./src/inference*.

##### Usage: 

```
java Main <options>
```

##### options:

- -data: the directory of data files
- -sigma: value of sigma; 1.0 by default
- -eps: value of eps; 0.8 by default
- -iter: value of maximum iterations for stochastic gradient ascent (in thousands); 5000 thousand by default
- -outer: value of maximum outer iterations for coordinate gradient ascent; 7 by default
- -verbose: whether see verbose output or not (0: show limited outputs; other: show all outputs); 1 by default
- -start: start timeslice; 120 by default
- -time: length of timespan; 100 by default

##### Example: 

```
  java Main -data ../../data/cosponsor/ -sigma 1 -eps 0.8 -iter 5000 -verbose 1 -start 120 -time 100
```


### Contact:
[Yupeng Gu](http://web.cs.ucla.edu/~ypgu/)

