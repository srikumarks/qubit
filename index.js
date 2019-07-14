import { complex, Complex } from './complex';
import { register } from './register';
import { Q0, Q1, QHP, QHM, qubit } from './qubit';
import {
    X, Y, Z, XX, YY, ZZ, H, QFT, oracle, controlled,
    measure, rot, phase, Rk, evolve, connect, shor
} from './operators';

let { C0, C1 } = Complex;

module.exports = {
    C0: C0,
    C1: C1,
    Q0: Q0,
    Q1: Q1,
    QHP: QHP,
    QHM: QHM,
    not: X,
    X: X,
    Y: Y,
    Z: Z,
    XX: XX,
    YY: YY,
    ZZ: ZZ,
    H: H,
    QFT: QFT,
    shor: shor,
    oracle: oracle,
    measure: measure,
    rot: rot,
    phase: phase,
    Rk: Rk,
    evolve: evolve,
    controlled: controlled,
    complex: complex,
    qubit: qubit,
    register: register,
    connect: connect,
    _: {
        r0: register({ctrl: Q0, data: Q1}),
        r1: register({ctrl: Q1, data: Q1}),
        bv: [register({x0: Q0, x1: Q0}), register({x0: Q1, x1: Q0}), register({x0: Q0, x1: Q1}), register({x0: Q1, x1: Q1})],
        bv4: (function () {
            let r = [];
            for (let i = 0; i < 16; ++i) {
                let b0 = (i & 1) ? Q1 : Q0;
                let b1 = ((i >> 1) & 1) ? Q1 : Q0;
                let b2 = ((i >> 2) & 1) ? Q1 : Q0;
                let b3 = ((i >> 3) & 1) ? Q1 : Q0;
                r.push(register({x0: b0, x1: b1, x2: b2, x3: b3}));
            }
            return r;
        }()),
        bv1: register({x0: Q1, x1: Q0, x2: Q1, x3: Q1}),
        ckt0: connect([
                    H('ctrl', 'data'),
                    controlled(X, 'ctrl', 'data')
                ]),
        deutsch: function (uf) {
            return connect([
                H('data'),
                uf('ctrl','data'),
                H('data')
            ]);   
        },
        qft2: QFT(['x0', 'x1']),
        qft4: QFT(['x0', 'x1', 'x2','x3'])
    }
}