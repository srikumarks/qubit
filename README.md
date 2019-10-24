# qubit

A silly little javascript based quantum computer simulator that I built as I
worked through the course https://coursera.org/learn/quantum-computing-algorithms

The course was fun to go through and this module has many operators you can
use and compose - such as n-bit Hadamard, QFT, etc. I haven't yet implemented
Shor's and Grover's algorithms. Coming soon. 

I saw some simulator libraries before this (not an exhaustive search), and I
realized I didn't want to index the qubits using integer indices. So this
library addresses the qubits using names. One consequence of this choice is
that the circuit can be defined independent of the particular bit order of
the representation chosen and will work exactly the same. This is in the
spirit of quantum mechanics, since QM doesn't care which order you want to
pack your system state into - i.e. which order you do the state kronecker
product. The physics remains the same. So I thought the code implementing
such a circuit must also show that representation independence.

The register you use to run your final circuit on can have any number of
bits. As long as it has qubits addressable using the names referenced in the
circuit, you're good to go. Of course, you must also provide inputs according
to that, so the construction of a register also permits you to use a "name ->
qubit" map instead of an array of qubits.

> **Note**: Needs testing! Since I was building this more to understand
> than as production code, I was fiddling around with representations 
> all the time without a priori design ... since I didn't have a clear
> idea of the design. Now I do have some idea and might rewrite and 
> expand this. I think it is more instructive to build quantum circuits
> in code than visually ... once the circuits start getting complex.

## How to use

```
# cd to qubit directory
node
> let Q = require('./start')
> let r = Q.register({x: Q.Q0, c: Q.Q1});
> r.toString()
'|0>|1>'
> let cnot = Q.controlled(Q.X, 'c', 'x'); // Qubit 'x' is controlled by qubit 'c'.
> let result = cnot(r);
> result.toString();
'|1>|1>'
> let cnotH = Q.connect([Q.H('c'), cnot]);
> result = cnotH(r);
> result.toString();
'(1/√2)|0⟩|0⟩ + (-1/√2)|1⟩|1⟩'
> result.expand().toString(); // State basis vector coefficients are (mag,phase_angle_degrees)
'(1/√2,0)|0> + (0,0)|1> + (0,0)|2> + (1/√2,180)|3>'
> r = Q.register({x: Q.Q1, c: Q.Q1});
> r.toString()
'|1>|1>'
> result = cnotH(r);
> result.toString();
'(1/√2)|1⟩|0⟩ + (-1/√2)|0⟩|1⟩'
```




