# COLA: a determinization, complementation, and containment checking library for Büchi automata

COLA has been built on top of [SPOT](https://spot.lrde.epita.fr/) and inspired by [Seminator](https://github.com/mklokocka/seminator).


### List of algorithms
Complementing UNBA (Under construction): the [slice-based algorithm](https://arxiv.org/abs/2005.09125v2)

Divide and conquer approach for determinization and complementaton of nondeterministic Büchi automata

Containment checking based on the determinization and complementation above

### Requirements
* [Spot](https://spot.lrde.epita.fr/)

```
./configure --enable-max-accsets=128
```
One can set the maximal number of colors for an automaton when configuring Spot with --enable-max-accsets=INT
```
make && make install
```

### Compilation
Please run the following steps to compile COLA after cloning this repo:
```
autoreconf -i
```
```
./configure [--with-spot=where/spot/is/installed]
```
```
make
```

Then you will get an executable file named **cola** !

### Determinization
Input an NBA from "filename", and run ```./cola --algo=cola filename```, then you will get an equivalent deterministic automaton on standard output!

To output the result to a file, use ```./cola --algo=cola filename -o out_filename```

To output a deterministic Parity automaton, use ```./cola --algo=cola filename```

To output a deterministic Rabin automaton, use ```./cola --algo=cola filename --rabin```

To output a complement automaton, use ```./cola --algo=cola filename --complement```

### Complementation
Complementation using decomposition: ```./cola --algo=comp [options] filename```

Options:
* ```--merge-iwa``` - merge all IWA SCCs
* ```--merge-det``` - merge all deterministic SCCs
* ```--verbose=1``` - outputs state names for debugging
* ```--tgba``` - outputs a TGBA with two colours
* ```--iw-sim``` - simulation on IWA SCCs
* ```--scc-compl``` - complementation for each SCC separately
