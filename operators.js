import { flatten } from './utils';
import { complex, Complex } from './complex';
import { register } from './register';
import { Q0, Q1, qubit } from './qubit';

let {C0, C1, CI, CIC} = Complex;

///////////////////////////////
// Basic operators.

function Gate(r_fn) {
    return function (...qs) {
        return apply(r_fn, flatten(qs));
    };
}

export let X = Gate(r_X);
export let Y = Gate(r_Y);
export let Z = Gate(r_Z);
export let H = Gate(r_Hn);
export let swap = Gate(r_swap);

function PGate(fn) {
    return function f(param, ...qs) {
        if (qs.length > 0) { return apply(fn(param), flatten(qs)); }
        return function (qs) {
            return f(param, ...qs);
        };
    };
}

export let XX = PGate(r_XX);
export let YY = PGate(r_YY);
export let ZZ = PGate(r_ZZ);
export let rot = PGate(r_rot);
export let phase = PGate(r_phase);
export let Rk = PGate(r_Rk);
export let evolve = PGate(r_evolve);

export function controlled(fn, c, ...qs) {
    if (c !== undefined) {
        if (qs.length > 0) {
            if (qs.indexOf(c) >= 0) {
                throw new Error('Control qubit has to be different from op qubit.');
            }
            let qfn = fn(qs);
            return function (reg) {
                if (qs.map(reg.lookup.bind(reg)).indexOf(reg.lookup(c)) >= 0) {
                    // Above we do the easy name checks. Here we do the
                    // index check, which requires the register name mapper.
                    throw new Error('Control qubit has to be different from op qubit.');
                }
                // If the control qbit is |0>, then it will leave
                // the entire register untouched. If the control
                // qbit is |1>, then it will apply the unitary transform
                // to the rest of the state. The net result is just a
                // superposition of these two conditions which we get
                // by concatenating the two partial result states.
                let reg0 = reg.project(c, 0);
                let qfnr = qfn(reg);
                let qfnr1 = qfnr.project(c, 1);
                return register(reg0.s.concat(qfnr1.s)).normalize().simplify().bindFrom(reg);
            };
        }
        return function (qs) {
            return controlled(fn, c, ...qs);
        };
    }
    return function (c_and_qs) {
        return controlled(fn, ...c_and_qs);
    };
}

// Quantum Fourier Transform
// We use the Rk decomposition here, except that at the end
// of the gate sequence, we swap the opposite qubits to retain
// the same qubit order for the register.
// xs is expected to be a qubit index array of the form
//  ['x0','x1', ... , 'xn']
// or the indices directly, given in "least significant qubit first"
// order.
export function QFT(xs) {
    let ops = [];
    for (let i = 0; i < xs.length; ++i) {
        for (let j = 1; j <= i; ++j) {
            ops.push(controlled(Rk(j+1), xs[i-j], xs[i]));
        }

        ops.push(H(xs[i]));
    }

    ops.reverse();

    for (let i = 0; i < xs.length; ++i) {
        let j = xs.length - i - 1;
        if (i >= j) { break; }
        ops.push(swap(xs[i], xs[j]));
    }

    return connect(ops);
}

// f :: Int -> Int     The function to turn into a quantum processor Q|x>|y> = |x>|y xor f(x)>
// xs :: [String]      The list of qubits (least first order) for x.
// ys :: [String]      The list of qubits (least first order) for y.
export function oracle(f, xs, ys) {
    let xys = {};
    for (let i = 0; i < xs.length; ++i) {
        xys[xs[i]] = true;
    }
    for (let i = 0; i < ys.length; ++i) {
        xys[ys[i]] = true;
    }
    if (Object.keys(xys).length < xs.length + ys.length) {
        throw new Error('xs and ys cannot have common qubits.');
    }

    return function (reg) {
        reg = reg.separate(...xs); // Guarantees that the xs qubits are all either Q0 or Q1 in each state component.

        let ixs = xs.map((x) => reg.lookup(x));
        let iys = ys.map((y) => reg.lookup(y));
        let nx = ixs.length, ny = iys.length;
        let Nx = 1 << nx, Ny = 1 << ny;

        for (let i = 0; i < this.s.length; ++i) {
            let sk = this.s[i].k;
            let x = 0;
            for (let ix = 0; ix < nx; ++ix) {
                let i = ixs[ix];
                if (sk[i] === Q1) {
                    x |= (1 << ix);
                }
            }
            let yBitsToFlip = f(x);
            for (let iy = 0; iy < ny; ++iy) {
                let i = iys[iy];
                let by = (yBitsToFlip >> iy) & 1;
                if (by) {
                    // If this should be a bit flip, then
                    // we exchange the zero and one amplitudes
                    // of the identified qubit.
                    let q = sk[i];
                    sk[i] = q === Q0 ? Q1 : (q === Q1 ? Q0 : qubit(q.one, q.zero));
                }
            }
        }
        return reg;
    };
}

// Puts the given qubits in an eigenstate of theirs
// depending on their amplitudes, at random.
export function measure(...xs) {
    return function (reg) {
        return reg.measure(...xs);
    };
}

// Prepares a "program" on a register that applies
// a sequence of operators to the register in the given
// order. An "operator" is a function that takes a register
// and produces a register as the result. The output could
// be a modified version of the input, so be prepared.
//
// The resultant circuit itself can be used as such an opertaor.
export function connect(opseq) {
    return function (reg) {
        for (let i = 0; i < opseq.length; ++i) {
            reg = opseq[i](reg);
        }
        return reg;
    };
}

export function shor(xs,f) {
    let ys = xs.map((x) => 'y_'+x);
    return connect([
        H(...xs),
        oracle(f,xs,ys),
        QFT(xs)
    ]);
}

////////////////////////////////
// Basic operator implementations

// Evolves the phase of only the |1> part of the qubit.
function r_phase(angle) {
    return function ([q]) {
        let phase = Complex.phase(angle);
        return register([{ c: C1, k: [qubit(q.zero, q.one.times(phase))] }]);
    };
}

function r_Rk(k) {
    return r_phase(2 * Math.PI * Math.pow(2, -k));
}

/*
function amplitude(qs, k) {
    let a = C1;
    for (let i = 0; i < qs.length; ++i) {
        a = a.times(k % 2 ? qs[i].one : qs[i].zero);
        k = Math.floor(k / 2);
    }
    return a;
}

function basisVec(k, n) {
    let s = [];
    for (let i = 0; i < n; ++i) {
        s.push(k % 2 ? Q1 : Q0);
        k = Math.floor(k/2);
    }
    return s;
}
*/



// Changes the complex phase of the state vector.
// This doesn't affect probabilities, but just changes
// the amplitude. Could be useful during superposition.
function r_evolve(angle) {
    return function ([q]) {
        let phase = Complex.phase(angle);
        return register({ c: C1, k: [qubit(q.zero.times(phase), q.one.times(phase))] });
    };
}

// Rotates a qubit in its state space.
function r_rot(angle) {
    return function ([q]) {
        let c = Math.cos(angle), s = Math.sin(angle);
        return register([{
            c: C1, k: [
                qubit(
                    q.zero.scale(c).plus(q.one.scale(-s)),
                    q.zero.scale(s).plus(q.one.scale(c))
                )
            ]
        }]);
    };
}

// Inverts a qubit.
function r_X([q]) {
    return register([{c:C1, k:[qubit(q.one,q.zero)]}]);
}

function r_Y([q]) {
    return register([{c:C1, k:[qubit(CIC.times(q.one), CI.times(q.zero))]}]);
}

function r_Z([q]) {
    return register([{c:C1, k:[qubit(q.zero, q.one.neg())]}]);
}

// Hadamard of multiple independent qubits. A common first
// step in quantum algorithms.
function r_Hn(qs) {
    return register([
        {c:C1, k: qs.map((q) => qubit(q.zero.plus(q.one), q.zero.minus(q.one)))}
    ]);
}

function r_XX(angle) {
    let phase = Complex.phase(angle).rot240();
    return function ([a,b]) {
        return register([
            {c: a.zero.times(b.zero).plus(phase.times(a.one).times(b.one)),
             k: [Q0,Q0]},
            {c: a.zero.times(b.one).plus(a.one.times(b.zero).rot240()),
             k: [Q0,Q1]},
            {c: a.zero.times(b.one).rot240().plus(a.one.times(b.zero)),
             k: [Q1,Q0]},
            {c: a.zero.times(b.zero).times(phase).plus(a.one.times(b.one)),
             k: [Q1,Q1]}
        ]);
    };
}

function r_YY(angle) {
    let c = complex(Math.cos(angle));
    let s = complex(0,Math.sin(angle));
    return function ([a,b]) {
        return register([
            {c: a.zero.times(b.zero).times(c).plus(s.times(a.one).times(b.one)),
             k: [Q0, Q0]},
            {c: a.zero.times(b.one).times(c).plus(s.neg().times(a.one.times(b.zero))),
             k: [Q0, Q1]},
            {c: a.zero.times(b.one).times(s.neg()).plus(c.times(a.one).times(b.zero)),
             k: [Q1, Q0]},
            {c: a.zero.times(b.zero).times(s).plus(c.times(a.one).times(b.one)),
             k: [Q1, Q1]}
        ]);
    };
}

function r_ZZ(angle) {
    let phase = Complex.phase(angle/2);
    let phasec = phase.conj();
    return function ([a,b]) {
        return register([
            {c: a.zero.times(b.zero).times(phase),
             k: [Q0,Q0]},
            {c: a.zero.times(b.one).times(phasec),
             k: [Q0,Q1]},
            {c: a.one.times(b.zero).times(phasec),
             k: [Q1,Q0]},
            {c: a.one.times(b.one).times(phase),
             k: [Q1,Q1]}
        ]);
    };
}

function r_swap([a,b]) {
    return register([{c:C1, k:[b,a]}]);
}

// Utility function to take an operator implementation and
// apply it to the given qubits of the register. The qubits
// identified by name ideally, but if not they need to be
// identified by index. So ixs is an array of names or indices
// (mixture is ok).
function apply(fn, ixs) {
    // Selects the qubits identified by ixs,
    // from the ks qubits array. The reg is used
    // for the name->index mapping.
    let sel = function (reg, ks) {
        let result = [];
        for (let i = 0; i < ixs.length; ++i) {
            let ix = (ixs[i] in reg.names ? reg.names[ixs[i]] : ixs[i]);
            result.push(ks[ix]);
        }
        return result;
    };

    // Applying an operator leaves the kronecker product sequence
    // unaffected in the register. This is an assumption needed for
    // the operator types. If the operator results in a one-to-one
    // change of qubit, then it is easy. If the operator is an
    // entangling operator, we need to split out the register's
    // states using the output of the operator so that each entry
    // in the result register uses base states. Hence the nested
    // loop. 
    //
    // Note that any of this isn't a viable simulation for even
    // a quantum system of more than tens of qubits. It would
    // require phenomenal amounts of RAM to work it. That's sort
    // of the whole point of and interest in quantum systems.
    return function (reg) {
        let s = [];
        for (let i = 0; i < reg.s.length; ++i) {
            let si = reg.s[i];
            let qs = fn(sel(reg, si.k));
            let part = [];
            for (let j = 0; j < qs.s.length; ++j) {
                let copy = [].concat(si.k);
                for (let m = 0; m < qs.s[j].k.length; ++m) {
                    let ix = (ixs[i] in reg.names ? reg.names[ixs[m]] : ixs[m]);
                    copy[ix] = qs.s[j].k[m];
                }
                part.push({
                    c: si.c.times(qs.s[j].c),
                    k: copy
                });
            }
            s.push(part);
        }
        return register(flatten(s)).bindFrom(reg);
    };
}
