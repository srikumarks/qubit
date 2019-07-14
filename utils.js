
export function showNum(n) {
    if (Math.abs(n - sq2inv) < eps) { return '1/√2'; }
    if (Math.abs(n + sq2inv) < eps) { return '-1/√2'; }
    if (Math.abs(Math.round(n) - n) < eps) { return ''+Math.round(n); }
    if (Math.abs(1/n - Math.round(1/n)) < eps) { return n > 0 ? ('1/'+Math.round(1/n)) : ('-1/'+Math.round(-1/n)); }
    n = Math.round(n * 10000) / 10000;
    return ''+n;
}

export function flatten(arrOfArr) {
    return Array.prototype.concat.apply([], arrOfArr);
}

export const sq2inv = 1.0 / Math.sqrt(2);
export const eps = 1e-6;
