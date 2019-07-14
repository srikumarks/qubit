// A simple library for creating and composing
// unitary transforms. Intention is to have it
// be easy to map between wiring diagrams and
// code so that eventually I can build a visual
// simulator.
//
// Why am I writing this? There are plenty of such
// code around, but I want to take this opportunity
// to think through for myself. That's basically it.
// As usual, I'm never satisfied until I implement
// shit myself.

import { eps } from './utils';
import { Complex } from './complex';

let { C0, C1 } = Complex;

function ket(state) {
    return '|'+state+'⟩';
}

function bra(state) {
    return '⟨'+state+'|';
}

function braket(l,r) {
    return '⟨'+l+'|'+r+'⟩';
}

export class Qubit {

    constructor(zero, one) {
        this.zero = zero;
        this.one = one;
    }

    toString(role) {
        let q = knownType(this);
        if (q) {
            return q;
        }
    
        if (role === 'coeff') {
            return '(' + showQubit(this) + ')';
        }
    
        return showQubit(this);
    }
}

function knownType(qubit) {
    if (qubit.zero.dist(Q0.zero) < eps && qubit.one.dist(Q0.one) < eps) {
        return Q0.form;
    }
    if (qubit.zero.dist(Q1.zero) < eps && qubit.one.dist(Q1.one) < eps) {
        return Q1.form;
    }
    if (qubit.zero.dist(QHP.zero) < eps && qubit.one.dist(QHP.one) < eps) {
        return QHP.form;
    }
    if (qubit.zero.dist(QHM.zero) < eps && qubit.one.dist(QHM.one) < eps) {
        return QHM.form;
    }
    return null;
}

export function qubit(zero, one) {
    let norm = 1.0 / Math.sqrt(zero.magsq() + one.magsq());
    return new Qubit(zero.scale(norm), one.scale(norm));
}

export let Q0 = qubit(C1, C0);
export let Q1 = qubit(C0, C1);
export let QHP = qubit(C1, C1);
export let QHM = qubit(C1, C1.neg());

Q0.form = ket(0);
Q1.form = ket(1);
QHP.form = ket('+');
QHM.form = ket('-');

function formAsString() { return this.form; }

Q0.toString = formAsString;
Q1.toString = formAsString;
QHP.toString = formAsString;
QHM.toString = formAsString;


function showQubit(q) {
    let s = '';
    if (q.zero.dist(C0) > eps) { s += q.zero.toString('coeff') + '|0>'; }
    if (q.one.dist(C0) > eps) { 
        if (s.length > 0) {
            s += ' + ';
        }
        s += q.one.toString('coeff') + '|1>';
    }
    return s;
}

