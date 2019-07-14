import { showNum, eps } from './utils';

export class Complex {
    constructor(re, im) {
        this.re = re;
        this.im = im;
    }

    dist(c) {
        let dx = this.re - c.re;
        let dy = this.im - c.im;
        return Math.sqrt(dx * dx + dy * dy);
    }

    conj() {
        return new Complex(this.re, -this.im);
    }

    rot90() {
        return new Complex(-this.im, this.re);
    }

    rot180() {
        return new Complex(-this.re, -this.im);
    }

    rot240() {
        return new Complex(this.im, -this.re);
    }

    plus(c) {
        return new Complex(this.re + c.re, this.im + c.im);   
    }

    scale(s) {
        return new Complex(this.re * s, this.im * s);
    }

    neg() {
        return this.rot180();
    }

    minus(c) {
        return new Complex(this.re - c.re, this.im - c.im);
    }

    times(c) {
        return new Complex(
            this.re * c.re - this.im * c.im,
            this.re * c.im + this.im * c.re
            );
    }

    divby(c) {
        return this.times(c.inv());
    }

    magsq() {
        return this.re * this.re + this.im * this.im;
    }

    mag() { 
        return Math.sqrt(this.magsq());
    }

    isSignificant() {
        return this.mag() > eps;
    }

    unit() {
        return this.scale(1.0 / this.mag());
    }

    inv() {
        return this.conj().scale(1.0 / this.magsq());
    }

    toString(role) {
        if (this.dist(Complex.C0) < eps) { return Complex.C0.toString(role); }
        if (this.dist(Complex.C1) < eps) { return Complex.C1.toString(role); }
        if (this.dist(Complex.CI) < eps) { return Complex.CI.toString(role); }
        if (this.dist(Complex.CIC) < eps) { return Complex.CIC.toString(role); }
        if (role === 'coeff') {
            return '(' + this.toString() + ')';
        }
        let resig = Math.abs(this.re) > eps;
        let imsig = Math.abs(this.im) > eps;
        return (resig ? showNum(this.re) : '') + (imsig ? (resig ? ' + ' : '') + showNum(this.im) + 'i' : '');
    } 

    static phase(angle) {
        return new Complex(Math.cos(angle), Math.sin(angle));
    }
}

Complex.C0 = new Complex(0,0);
Complex.C1 = new Complex(1,0);
Complex.CI = new Complex(0,1);
Complex.CIC = new Complex(0,-1);
Complex.zero = Complex.C0;
Complex.unity = Complex.C1;
Complex.one = Complex.C1;
Complex.i = Complex.CI;

Complex.C0.toString = function () { return '0'; };
Complex.C1.toString = function (role) { return role === 'coeff' ? '' : '1'; };
Complex.CI.toString = function () { return 'i'; };
Complex.CIC.toString = function () { return '-i'; };

export function complex(re, im) {
    if (re && re.conj) { return re; }
    return new Complex(re, im);
}

complex.phase = Complex.phase;
