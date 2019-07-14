import { showNum, eps } from './utils';
import { Complex } from './complex';
import { Q0, Q1 } from './qubit';
let { C1, C0 } = Complex;

export class Register {

    constructor(state) {
        this.n = state.length > 0 ? state[0].k.length : 0;
        this.s = state;
        this.names = {};
    }

    // Binds a name to a given register index, or other name.
    name(name, ref) {
        this.names[name] = this.lookup(ref);
        return this;
    }

    lookup(ref) {
        return (ref in this.names) ? this.names[ref] : ref;
    }

    lookupAll(refs) {
        return refs.map((r) => r in this.names ? this.names[r] : r);
    }

    // Binds given names to register qubits. The nameMap is
    // an object mapping new names to existing names or indices.
    bind(nameMap) {
        for (let k in nameMap) {
            this.names[k] = this.lookup(nameMap[k]);
        }
        return this;
    }

    bindFrom(reg) {
        if (this.n !== reg.n) {
            throw new Error('Incompatible registers to transfer names.');
        }
        this.bind(reg.names);
        return this;
    }

    // Some opportunity for tree structuring and lazy
    // expansion here.
    kron(reg) {
        let s = [];
        for (let i = 0; i < this.s.length; ++i) {
            let si = this.s[i];
            for (let j = 0; j < reg.s.length; ++j) {
                let sj = reg.s[j];
                s.push({ c: si.c.times(sj.c), k: si.k.concat(sj.k) });
            }
        }
        let r = new Register(s);
        for (let k in reg.names) {
            if (k in this.names) {
                throw new Error('Clash in qubit names. Names need to be unique.');
            }
            r.names[k] = this.n + reg.names[k];
        }
        return r;
    }

    superpose(cself, creg, reg) {
        if (this.n !== reg.n) {
            throw new Error('Number of states must be the same for superposition.');
        }

        // Check name compatibility. If same name is present
        // in both, they must refer to the same qubit.
        for (let n in reg.names) {
            if (n in this.names) {
                if (this.names[n] !== reg.names[n]) {
                    throw new Error('Name collision: Name "' + n + '" means different things in superposed registers.');
                }
            }
        }
        for (let n in this.names) {
            if (n in reg.names) {
                if (this.names[n] !== reg.names[n]) {
                    throw new Error('Name collision: Name "' + n + '" means different things in superposed registers.');
                }
            }
        }

        let scaler = (c) => (e) => { return { c: e.c.times(c), k: e.k }; };

        let s = this.s.map(scaler(cself)).concat(reg.s.map(scaler(creg)));
        return register(s).bindFrom(this).bindFrom(reg); // Use either name as long as no collision.
    }

    // Splits each state vector component so that the qubits
    // identified by the xs are all either Q0 or Q1 in each
    // component.
    separate(...xs) {
        let ixs = this.lookupAll(xs);
        for (let i = 0; i < ixs.length; ++i) {
            let xi = ixs[i];
            let states = [];
            for (let si = 0; si < this.s.length; ++si) {
                let s = this.s[si];
                if (s.k[xi].zero.mag() > eps) {
                    states.push({ c: s.c.times(s.k[xi].zero), k: s.k.map((q, qi) => qi === xi ? Q0 : q) });
                }
                if (s.k[xi].one.mag() > eps) {
                    states.push({ c: s.c.times(s.k[xi].one), k: s.k.map((q, qi) => qi === xi ? Q1 : q) });
                }
            }
            this.s = states;
        }
        return this;
    }

    // When the indicated qubits are measured, they fall into
    // one of their eigenstates depending on their amplitudes,
    // at random.
    measure(...xs) {
        let ixs = this.lookupAll(xs);

        for (let i = 0; i < ixs.length; ++i) {
            let ix = ixs[i];
            let amp0 = C0, amp1 = C0;
            for (let si = 0; si < this.s.length; ++si) {
                let s = this.s[si];
                amp0 = amp0.plus(s.c.times(s.k[ix].zero));
                amp1 = amp1.plus(s.c.times(s.k[ix].one));
            }
            let p0 = amp0.magsq();
            let p1 = amp1.magsq();
            if (Math.abs(p0 + p1 - 1) > eps) {
                throw new Error('BUG: Some calculation error here.');
            }

            let choice = Math.random() <= p0 ? Q0 : Q1;
            for (let si = 0; si < this.s.length; ++si) {
                this.s[si].k[ix] = choice;
            }
        }

        return this;
    }


    // name identifies a qubit and state is either 0 or 1.
    // Returns an un-normalized register so that you can
    // calculate its probability before you normalize it.
    project(name, state) {
        let ix = this.lookup(name);
        let s = this.s.map((e) => {
            let qix = e.k[ix];
            return { c: e.c.times(state ? qix.one : qix.zero), k: e.k.map((q, i) => i === ix ? (state ? Q1 : Q0) : q) };
        });
        return (new Register(s)).bindFrom(this);
    }

    // Given state vector calculates its sum of squared amplitudes.
    // This is useful for projections and normalization.
    probability() {
        let c = 0;
        for (let i = 0; i < this.s.length; ++i) {
            c += stateVec[i].c.magsq();
        }
        return c;
    }

    // Normalizes the state amplitudes of this register.
    // A register created using register() is already
    // normalized, but on occasion (ex: reg.project()), 
    // you get an unnormalized register.
    normalize() {
        let n = 0, s = this.s;
        for (let i = 0; i < s.length; ++i) {
            n += s[i].c.magsq();
        }

        if (Math.sqrt(n) < eps) {
            throw new Error('Invalid register with near zero probability.');
        }

        n = 1.0 / n;
        for (let i = 0; i < s.length; ++i) {
            s[i].c = s[i].c.scale(n);
        }

        return this;
    }

    simplify() {
        let fstate = [];
        for (let i = 0; i < this.s.length; ++i) {
            if (this.s[i].c.mag() > eps) {
                fstate.push(this.s[i]);
            }
        }
        this.s = fstate;
        return this;
    }

    expand() {
        if (this.n > 24) {
            throw new Error('Too big quantum state for this poor classical machine.')
        }

        let amplitudes = new Float64Array(2 << this.n);

        let N = 1 << this.n;
        for (let i = 0; i < this.s.length; ++i) {
            let s = this.s[i];
            let k = s.k;
            for (let j = 0; j < N; ++j) {
                let { re, im } = s.c;
                for (let m = 0; m < k.length; ++m) {
                    let q = ((j >> m) & 1) ? k[m].one : k[m].zero;
                    let _re = re * q.re - im * q.im;
                    let _im = re * q.im + q.re * im;
                    re = _re;
                    im = _im;
                }
                amplitudes[j * 2] += re;
                amplitudes[j * 2 + 1] += im;
            }
        }

        amplitudes.toString = showExpandedAmplitudes;

        return amplitudes;
    }

    toString() {
        let s = '';
        let count = 0;
        for (let i = 0; i < this.s.length; ++i) {
            if (count > 0) { s += ' + ' };
            if (this.s[i].c.dist(C0) > eps) {
                let frag = this.s[i].c.toString('coeff') + this.s[i].k.map((q) => q.toString('coeff')).join('');
                s += frag;
                ++count;
            }
        }
        return s;
    }

}


function showExpandedAmplitudes() {
    let amplitudes = this;
    let s = [];
    for (let i = 0; i < amplitudes.length; i += 2) {
        let magsq = amplitudes[i] * amplitudes[i] + amplitudes[i+1] * amplitudes[i+1];
        let phase = Math.round(Math.atan2(amplitudes[i+1], amplitudes[i]) * 180 / Math.PI);
        s.push('(' + showNum(Math.sqrt(magsq)) + ',' + showNum(phase) + ')|' + (i >> 1) + '>');
    }
    return s.join(' + ');
}


// state is an array of {c: complex, k: kron}.
export function register(state) {
    if (state instanceof Array) {
        return (new Register(state)).normalize();
    }
    return register.eigen(state);
}

register.eigen = function (qs) {
    if (qs instanceof Array) {
        return register([{c:C1,k:qs}]);
    }
    let ks = [], names = {}, i = 0;
    for (let name in qs) {
        ks.push(qs[name]);
        names[name] = i;
        ++i;
    }
    return register([{c:C1,k:ks}]).bind(names);
};
