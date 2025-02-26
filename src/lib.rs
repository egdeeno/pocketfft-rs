#![allow(non_camel_case_types)]
use std::ptr::null_mut;
use std::slice::from_raw_parts_mut;

const NFCT: usize = 25;

#[repr(C)]
pub struct cfft_plan_i {
    pub packplan: cfftp_plan,
    pub blueplan: fftblue_plan,
}
pub type cfft_plan = *mut cfft_plan_i;

#[repr(C)]
pub struct fftblue_plan_i {
    pub n: usize,
    pub n2: usize,
    pub plan: cfftp_plan,
    pub mem: Vec<f64>,
}
pub type fftblue_plan = *mut fftblue_plan_i;

#[repr(C)]
pub struct cfftp_plan_i {
    pub length: usize,
    pub nfct: usize,
    pub mem: Vec<cmplx>,
    pub fct: [cfftp_fctdata; NFCT],
}
pub type cfftp_plan = *mut cfftp_plan_i;

#[derive(Copy, Clone)]
#[repr(C)]
pub struct cfftp_fctdata {
    pub fct: usize,
    pub tw_index: usize,
    pub tws_index: usize,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct cmplx {
    pub r: f64,
    pub i: f64,
}
impl PartialEq for cmplx {
    fn eq(&self, other: &Self) -> bool {
        self.r == other.r && self.i == other.i
    }
}

#[repr(C)]
pub struct rfft_plan_i {
    pub packplan: rfftp_plan,
    pub blueplan: fftblue_plan,
}
pub type rfft_plan = *mut rfft_plan_i;
#[repr(C)]
pub struct rfftp_plan_i {
    pub length: usize,
    pub nfct: usize,
    pub mem: Vec<f64>,
    pub fct: [rfftp_fctdata; NFCT],
}
pub type rfftp_plan = *mut rfftp_plan_i;

#[derive(Copy, Clone)]
#[repr(C)]
pub struct rfftp_fctdata {
    pub fct: usize,
    pub tw_index: usize,
    pub tws_index: usize,
}
//====================================
pub struct CfftPlan {
    pub packplan: CfftpPlan,
    pub blueplan: FftbluePlan,
}

pub struct RfftPlan {
    pub packplan: RfftpPlan,
    pub blueplan: FftbluePlan,
}

pub struct RfftpPlan {
    pub length: usize,
    pub nfct: usize,
    pub mem: Vec<f64>,
    pub fct: [RfftpFctdata; NFCT],
}
#[derive(Copy, Clone)]
pub struct RfftpFctdata {
    pub fct: usize,
    pub tw_index: usize,
    pub tws_index: usize,
}

#[derive(Default)]
pub struct FftbluePlan {
    pub n: usize,
    pub n2: usize,
    pub plan: CfftpPlan,
    pub mem: Vec<f64>,
}

#[derive(Default)]
pub struct CfftpPlan {
    pub length: usize,
    pub nfct: usize,
    pub mem: Vec<cmplx>,
    pub fct: [CfftpFctdata; NFCT],
}

#[derive(Default, Copy, Clone)]
pub struct CfftpFctdata {
    pub fct: usize,
    pub tw_index: usize,
    pub tws_index: usize,
}
//====================================

fn sincosm1pi(a: f64, res: &mut [f64]) {
    let mut s = a * a;
    // Approximate cos(pi*x)-1 for x in [-0.25,0.25]
    let mut r: f64 = -1.0369917389758117e-4;
    r = r.mul_add(s, 1.9294935641298806e-3);

    r = r.mul_add(s, -2.5806887942825395e-2);

    r = r.mul_add(s, 2.353_306_302_832_821e-1);

    r = r.mul_add(s, -1.335_262_768_853_800_6);

    r = r.mul_add(s, 4.058_712_126_416_762);

    r = r.mul_add(s, -4.934_802_200_544_679);

    let c = r * s;

    // Approximate sin(pi*x) for x in [-0.25,0.25]
    r = 4.6151442520157035e-4;
    r = r.mul_add(s, -7.3700183130883555e-3);

    r = r.mul_add(s, 8.214_586_894_932_394e-2);

    r = r.mul_add(s, -5.992_645_289_321_492e-1);

    r = r.mul_add(s, 2.550_164_039_873_269);

    r = r.mul_add(s, -5.167_712_780_049_952);

    r = r * s * a;
    s = a.mul_add(3.141_592_653_589_793, r);
    res[0] = c;
    res[1] = s;
}

fn calc_first_octant(den: usize, res: &mut [f64]) {
    let n = (den + 4) >> 3;
    if n == 0 {
        return;
    }
    res[0] = 1.0;
    res[1] = 0.0;
    if n != 1 {
        let l1 = (n as f64).sqrt() as usize;

        for i in 1..l1 {
            let a = (2.0 * i as f64) / (den as f64);
            sincosm1pi(a, &mut res[(2 * i)..]);
        }

        let mut start = l1;

        while start < n {
            let mut cs: [f64; 2] = [0.0; 2];
            let a = (2.0 * (start as f64)) / (den as f64);
            sincosm1pi(a, &mut cs);
            res[2 * start] = cs[0] + 1.0;
            res[2 * start + 1] = cs[1];
            let mut end = l1;
            if start + end > n {
                end = n - start;
            }
            for k in 1..end {
                let mut csx: [f64; 2] = [0.0; 2];
                csx[0] = res[2 * k];
                csx[1] = res[2 * k + 1];
                res[2 * (start + k)] = ((cs[0] * csx[0] - cs[1] * csx[1] + cs[0]) + csx[0]) + 1.0;
                res[2 * (start + k) + 1] = (cs[0] * csx[1] + cs[1] * csx[0]) + cs[1] + csx[1];
            }
            start += l1;
        }
        for i in 1..l1 {
            res[2 * i] += 1.0;
        }
    }
}

fn calc_first_quadrant(n: usize, res: &mut [f64]) {
    //assert_eq!(n, res.len());
    calc_first_octant(n << 1, &mut res[n..]);
    let ndone = (n + 2) >> 2;
    let mut i: usize = 0;
    let mut idx1: usize = 0;
    let mut idx2: usize = 2 * ndone - 2;
    while (i + 1) < ndone {
        res[idx1] = res[n + 2 * i];
        res[idx1 + 1] = res[n + 2 * i + 1];
        res[idx2] = res[n + 2 * i + 3];
        res[idx2 + 1] = res[n + 2 * i + 2];
        i += 2;
        idx1 += 2;
        idx2 -= 2;
    }
    if i != ndone {
        res[idx1] = res[n + 2 * i];
        res[idx1 + 1] = res[n + 2 * i + 1];
    }
}

fn calc_first_half(n: usize, res: &mut [f64]) {
    let ndone = (n + 1) >> 1;
    calc_first_octant(n << 2, &mut res[(n - 1)..]);
    let mut i4 = 0;
    let mut i = 0;
    while i4 < n && i4 <= (n - i4) {
        res[2 * i] = res[n - 1 + 2 * i4];
        res[2 * i + 1] = res[n - 1 + 2 * i4 + 1];
        i += 1;
        i4 += 4;
    }
    while (i4 as isize) - (n as isize) <= 0 {
        let xm = n - i4;
        res[2 * i] = res[n - 1 + 2 * xm + 1];
        res[2 * i + 1] = res[n - 1 + 2 * xm];
        i += 1;
        i4 += 4;
    }
    while i4 <= 3 * n - i4 {
        let xm = i4 - n;
        res[2 * i] = -1.0 * res[n - 1 + 2 * xm + 1];
        res[2 * i + 1] = res[n - 1 + 2 * xm];
        i += 1;
        i4 += 4;
    }
    while i < ndone {
        let xm = 2 * n - i4;
        res[2 * i] = -1.0 * res[n - 1 + 2 * xm];
        res[2 * i + 1] = res[n - 1 + 2 * xm + 1];
        i += 1;
        i4 += 4;
    }
}

fn fill_first_quadrant(n: usize, res: &mut [f64]) {
    let hsqt2 = 0.707_106_781_186_547_6;
    let quart = n >> 2;
    if (n & 7) == 0 {
        res[quart + 1] = hsqt2;
        res[quart] = hsqt2
    }
    let mut i = 2;
    let mut j = 2 * quart - 2;
    while i < quart {
        res[j] = res[i + 1];
        res[j + 1] = res[i];
        i += 2;
        j -= 2;
    }
}

fn fill_first_half(n: usize, res: &mut [f64]) {
    let half = n >> 1;
    if (n & 3) == 0 {
        let mut i = 0;
        while i < half {
            res[i + half] = -res[i + 1];
            res[i + half + 1] = res[i];
            i += 2;
        }
    } else {
        let mut i = 2;
        let mut j = 2 * half - 2;
        while i < half {
            res[j] = -res[i];
            res[j + 1] = res[i + 1];
            i += 2;
            j -= 2;
        }
    }
}

fn fill_second_half(n: usize, res: &mut [f64]) {
    if (n & 1) == 0 {
        for i in 0..n {
            res[i + n] = -res[i];
        }
    } else {
        let mut i = 2;
        let mut j = 2 * n - 2;
        while i < n {
            res[j] = res[i];
            res[j + 1] = -res[i + 1];
            i += 2;
            j -= 2;
        }
    }
}

fn sincos_2pibyn_half(n: usize, res: &mut [f64]) {
    if (n & 3) == 0 {
        calc_first_octant(n, res);
        fill_first_quadrant(n, res);
        fill_first_half(n, res);
    } else if (n & 1) == 0 {
        calc_first_quadrant(n, res);
        fill_first_half(n, res);
    } else {
        calc_first_half(n, res);
    }
}

fn sincos_2pibyn(n: usize, res: &mut [f64]) {
    sincos_2pibyn_half(n, res);
    fill_second_half(n, res);
}

fn largest_prime_factor(n: usize) -> usize {
    let mut n_temp = n;
    let mut res: usize = 1;
    let mut tmp = n_temp >> 1;
    while (tmp << 1) == n_temp {
        res = 2;
        n_temp = tmp;
        tmp = n_temp >> 1;
    }

    let mut limit = ((n_temp as f64) + 0.01).sqrt() as usize;
    let mut x = 3;
    while x <= limit {
        tmp = n_temp / x;
        while tmp * x == n_temp {
            res = x;
            n_temp = tmp;
            limit = ((n_temp as f64) + 0.01).sqrt() as usize;
        }
        x += 2;
    }
    if n_temp > 1 {
        res = n_temp;
    }

    return res;
}

fn cost_guess(n: usize) -> f64 {
    let lfp: f64 = 1.1; // penalty for non-hardcoded larger factors
    let mut n_temp = n;
    let ni = n;
    let mut result: f64 = 0.0;
    let mut tmp = n_temp >> 1;
    while (tmp << 1) == n_temp {
        result += 2.0;
        n_temp = tmp;
        tmp = n_temp >> 1;
    }

    let mut limit = ((n_temp as f64) + 0.01).sqrt() as usize;
    let mut x: usize = 3;
    while x <= limit {
        tmp = n_temp / x;
        while (tmp * x) == n_temp {
            if x <= 5 {
                result += x as f64;
            } else {
                result += lfp * (x as f64);
            }
            // penalize larger prime factors
            n_temp = tmp;
            limit = ((n_temp as f64) + 0.01).sqrt() as usize;
            tmp = n_temp / x;
        }
        x += 2;
    }
    if n_temp > 1 {
        if n_temp <= 5 {
            result += n_temp as f64;
        } else {
            result += lfp * (n_temp as f64);
        }
    }

    return result * (ni as f64);
}

// returns the smallest composite of 2, 3, 5, 7 and 11 which is >= n
fn good_size(n: usize) -> usize {
    if n <= 6 {
        return n;
    }

    let mut bestfac: usize = 2 * n;
    let mut f2: usize = 1;
    while f2 < bestfac {
        let mut f23 = f2;
        while f23 < bestfac {
            let mut f235 = f23;
            while f235 < bestfac {
                let mut f2357 = f235;
                while f2357 < bestfac {
                    let mut f235711 = f2357;
                    while f235711 < bestfac {
                        if f235711 >= n {
                            bestfac = f235711;
                        }
                        f235711 *= 11;
                    }
                    f2357 *= 7;
                }
                f235 *= 5;
            }
            f23 *= 3;
        }
        f2 *= 2;
    }
    return bestfac;
}

fn pass2b(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx]) {
    let cdim: usize = 2;

    if ido == 1 {
        for k in 0..l1 {
            ch[ido * k].r = cc[ido * cdim * k].r + cc[ido * (1 + cdim * k)].r;
            ch[ido * k].i = cc[ido * cdim * k].i + cc[ido * (1 + cdim * k)].i;
            ch[ido * (k + l1)].r = cc[ido * cdim * k].r - cc[ido * (1 + cdim * k)].r;
            ch[ido * (k + l1)].i = cc[ido * cdim * k].i - cc[ido * (1 + cdim * k)].i;
        }
    } else {
        for k in 0..l1 {
            ch[ido * k].r = cc[ido * cdim * k].r + cc[ido * (1 + cdim * k)].r;
            ch[ido * k].i = cc[ido * cdim * k].i + cc[ido * (1 + cdim * k)].i;

            ch[ido * (k + l1)].r = cc[ido * cdim * (k)].r - cc[ido * (1 + cdim * k)].r;
            ch[ido * (k + l1)].i = cc[ido * cdim * k].i - cc[ido * (1 + cdim * k)].i;
            for i in 1..ido {
                let mut t: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
                ch[i + ido * k].r = cc[i + ido * cdim * k].r + cc[i + ido * (1 + cdim * k)].r;
                ch[i + ido * k].i = cc[i + ido * cdim * k].i + cc[i + ido * (1 + cdim * k)].i;
                t.r = cc[i + ido * cdim * k].r - cc[i + ido * (1 + cdim * k)].r;
                t.i = cc[i + ido * cdim * k].i - cc[i + ido * (1 + cdim * k)].i;

                ch[i + ido * (k + l1)].r = wa[i - 1].r * t.r - wa[i - 1].i * t.i;
                ch[i + ido * (k + l1)].i = wa[i - 1].r * t.i + wa[i - 1].i * t.r;
            }
        }
    }
}

fn pass2f(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx]) {
    let cdim: usize = 2;

    if ido == 1 {
        for k in 0..l1 {
            ch[ido * k].r = cc[ido * cdim * k].r + cc[ido * (1 + cdim * k)].r;
            ch[ido * k].i = cc[ido * cdim * k].i + cc[ido * (1 + cdim * k)].i;
            ch[ido * (k + l1)].r = cc[ido * cdim * k].r - cc[ido * (1 + cdim * k)].r;
            ch[ido * (k + l1)].i = cc[ido * cdim * k].i - cc[ido * (1 + cdim * k)].i;
        }
    } else {
        for k in 0..l1 {
            ch[ido * k].r = cc[ido * cdim * k].r + cc[ido * (1 + cdim * k)].r;
            ch[ido * k].i = cc[ido * cdim * k].i + cc[ido * (1 + cdim * k)].i;
            ch[ido * (k + l1)].r = cc[ido * cdim * k].r - cc[ido * (1 + cdim * k)].r;
            ch[ido * (k + l1)].i = cc[ido * cdim * k].i - cc[ido * (1 + cdim * k)].i;
            for i in 1..ido {
                let mut t: cmplx = cmplx { r: 0.0, i: 0.0 };

                ch[i + ido * k].r = cc[i + ido * cdim * k].r + cc[i + ido * (1 + cdim * k)].r;
                ch[i + ido * k].i = cc[i + ido * cdim * k].i + cc[i + ido * (1 + cdim * k)].i;
                t.r = cc[i + ido * cdim * k].r - cc[i + ido * (1 + cdim * k)].r;
                t.i = cc[i + ido * cdim * k].i - cc[i + ido * (1 + cdim * k)].i;

                ch[i + ido * (k + l1)].r = wa[i - 1].r * t.r + wa[i - 1].i * t.i;
                ch[i + ido * (k + l1)].i = wa[i - 1].r * t.i - wa[i - 1].i * t.r;
            }
        }
    }
}

fn pass3b(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx]) {
    let cdim: usize = 3;
    let tw1r: f64 = -0.5;
    let tw1i: f64 = 0.866_025_403_784_438_6;

    if ido == 1 {
        for k in 0..l1 {
            //PREP3(0)
            let t0: cmplx = cc[ido * cdim * k]; //&cc[(0) + ido * ((0) + cdim * (k))];
            let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
            let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };

            t1.r = cc[ido * (1 + cdim * k)].r + cc[ido * (2 + cdim * k)].r;
            t1.i = cc[ido * (1 + cdim * k)].i + cc[ido * (2 + cdim * k)].i;
            t2.r = cc[ido * (1 + cdim * k)].r - cc[ido * (2 + cdim * k)].r;
            t2.i = cc[ido * (1 + cdim * k)].i - cc[ido * (2 + cdim * k)].i;
            ch[ido * k].r = t0.r + t1.r;
            ch[ido * k].i = t0.i + t1.i;

            let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
            let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
            ca.r = t0.r + tw1r * t1.r;
            ca.i = t0.i + tw1r * t1.i;
            cb.i = tw1i * t2.r;
            cb.r = -1.0 * (tw1i * t2.i);

            ch[ido * (k + l1)].r = ca.r + cb.r;
            ch[ido * (k + l1)].i = ca.i + cb.i;
            ch[ido * (k + l1 * 2)].r = ca.r - cb.r;
            ch[ido * (k + l1 * 2)].i = ca.i - cb.i;
        }
    } else {
        for k in 0..l1 {
            {
                //PREP3(0)
                let t0: cmplx = cc[ido * cdim * k];
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };

                //PMC (t1,t2,CC(idx,1,k),CC(idx,2,k))
                t1.r = cc[ido * (1 + cdim * k)].r + cc[ido * (2 + cdim * k)].r;
                t1.i = cc[ido * (1 + cdim * k)].i + cc[ido * (2 + cdim * k)].i;
                t2.r = cc[ido * (1 + cdim * k)].r - cc[ido * (2 + cdim * k)].r;
                t2.i = cc[ido * (1 + cdim * k)].i - cc[ido * (2 + cdim * k)].i;
                //CH(idx,k,0).r=t0.r+t1.r;
                ch[ido * k].r = t0.r + t1.r;

                //CH(idx,k,0).i=t0.i+t1.i;
                ch[ido * k].i = t0.i + t1.i;

                //PARTSTEP3a(1,2,tw1r,tw1i)
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
                ca.r = t0.r + tw1r * t1.r;
                ca.i = t0.i + tw1r * t1.i;
                cb.i = tw1i * t2.r;
                cb.r = -1.0 * (tw1i * t2.i);
                ch[ido * (k + l1)].r = ca.r + cb.r;
                ch[ido * (k + l1)].i = ca.i + cb.i;
                ch[ido * (k + l1 * 2)].r = ca.r - cb.r;
                ch[ido * (k + l1 * 2)].i = ca.i - cb.i;
            }
            for i in 1..ido {
                //PREP3(i)
                let t0: cmplx = cc[i + ido * cdim * k];
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 }; //cmplx { r: 0.0, i: 0.0 };
                t1.r = cc[i + ido * (1 + cdim * k)].r + cc[i + ido * (2 + cdim * k)].r;
                t1.i = cc[i + ido * (1 + cdim * k)].i + cc[i + ido * (2 + cdim * k)].i;
                t2.r = cc[i + ido * (1 + cdim * k)].r - cc[i + ido * (2 + cdim * k)].r;
                t2.i = cc[i + ido * (1 + cdim * k)].i - cc[i + ido * (2 + cdim * k)].i;
                ch[i + ido * k].r = t0.r + t1.r;
                ch[i + ido * k].i = t0.i + t1.i;

                //PARTSTEP3b(1,2,tw1r,tw1i)
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t0.r + tw1r * t1.r;
                ca.i = t0.i + tw1r * t1.i;
                cb.i = tw1i * t2.r;
                cb.r = -(tw1i * t2.i);
                da.r = ca.r + cb.r;
                da.i = ca.i + cb.i;
                db.r = ca.r - cb.r;
                db.i = ca.i - cb.i;

                ch[i + ido * (k + l1)].r = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.r
                    - wa[i - 1 + (1 - 1) * (ido - 1)].i * da.i;
                ch[i + ido * (k + l1)].i = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.i
                    + wa[i - 1 + (1 - 1) * (ido - 1)].i * da.r;

                ch[i + ido * (k + l1 * 2)].r =
                    wa[i - 1 + (ido - 1)].r * db.r - wa[i - 1 + (ido - 1)].i * db.i;
                ch[i + ido * (k + l1 * 2)].i =
                    wa[i - 1 + (ido - 1)].r * db.i + wa[i - 1 + (ido - 1)].i * db.r;
            }
        }
    }
}

fn pass3f(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx]) {
    let cdim: usize = 3;
    let tw1r: f64 = -0.5;
    let tw1i: f64 = -0.866_025_403_784_438_6;

    if ido == 1 {
        for k in 0..l1 {
            let t0: cmplx = cc[ido * cdim * k];
            let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
            {
                t1.r = cc[ido * (1 + cdim * k)].r + cc[ido * (2 + cdim * k)].r;
                t1.i = cc[ido * (1 + cdim * k)].i + cc[ido * (2 + cdim * k)].i;
                t2.r = cc[ido * (1 + cdim * k)].r - cc[ido * (2 + cdim * k)].r;
                t2.i = cc[ido * (1 + cdim * k)].i - cc[ido * (2 + cdim * k)].i;
            }
            ch[ido * k].r = t0.r + t1.r;
            ch[ido * k].i = t0.i + t1.i;
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t0.r + tw1r * t1.r;
                ca.i = t0.i + tw1r * t1.i;
                cb.i = tw1i * t2.r;
                cb.r = -(tw1i * t2.i);
                {
                    ch[ido * (k + l1)].r = ca.r + cb.r;
                    ch[ido * (k + l1)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 2)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 2)].i = ca.i - cb.i;
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let t0: cmplx = cc[ido * cdim * (k)];
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t1.r = cc[ido * (1 + cdim * k)].r + cc[ido * (2 + cdim * k)].r;
                    t1.i = cc[ido * (1 + cdim * k)].i + cc[ido * (2 + cdim * k)].i;
                    t2.r = cc[ido * (1 + cdim * k)].r - cc[ido * (2 + cdim * k)].r;
                    t2.i = cc[ido * (1 + cdim * k)].i - cc[ido * (2 + cdim * k)].i;
                }
                ch[ido * k].r = t0.r + t1.r;
                ch[ido * k].i = t0.i + t1.i;
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw1r * t1.r;
                    ca.i = t0.i + tw1r * t1.i;
                    cb.i = tw1i * t2.r;
                    cb.r = -(tw1i * t2.i);
                    {
                        ch[ido * (k + l1)].r = ca.r + cb.r;
                        ch[ido * (k + l1)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 2)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 2)].i = ca.i - cb.i;
                    }
                }
            }
            for i in 1..ido {
                let t0: cmplx = cc[(i) + ido * cdim * (k)];
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t1.r = cc[i + ido * (1 + cdim * k)].r + cc[i + ido * (2 + cdim * k)].r;
                    t1.i = cc[i + ido * (1 + cdim * k)].i + cc[i + ido * (2 + cdim * k)].i;
                    t2.r = cc[i + ido * (1 + cdim * k)].r - cc[i + ido * (2 + cdim * k)].r;
                    t2.i = cc[i + ido * (1 + cdim * k)].i - cc[i + ido * (2 + cdim * k)].i;
                }
                ch[(i) + ido * k].r = t0.r + t1.r;
                ch[(i) + ido * k].i = t0.i + t1.i;
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw1r * t1.r;
                    ca.i = t0.i + tw1r * t1.i;
                    cb.i = tw1i * t2.r;
                    cb.r = -(tw1i * t2.i);
                    {
                        da.r = ca.r + cb.r;
                        da.i = ca.i + cb.i;
                        db.r = ca.r - cb.r;
                        db.i = ca.i - cb.i;
                    }
                    {
                        ch[i + ido * (k + l1)].r = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.r
                            + wa[i - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1)].i = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.i
                            - wa[i - 1 + (1 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 2)].r =
                            wa[i - 1 + (ido - 1)].r * db.r + wa[i - 1 + (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 2)].i =
                            wa[i - 1 + (ido - 1)].r * db.i - wa[i - 1 + (ido - 1)].i * db.r;
                    }
                }
            }
        }
    }
}

fn pass4b(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx]) {
    let cdim: usize = 4;

    if ido == 1 {
        for k in 0..l1 {
            let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
            {
                t2.r = cc[ido * cdim * k].r + cc[ido * (2 + cdim * k)].r;
                t2.i = cc[ido * cdim * k].i + cc[ido * (2 + cdim * k)].i;
                t1.r = cc[ido * cdim * k].r - cc[ido * (2 + cdim * k)].r;
                t1.i = cc[ido * cdim * k].i - cc[ido * (2 + cdim * k)].i;
            }
            {
                t3.r = cc[ido * (1 + cdim * k)].r + cc[ido * (3 + cdim * k)].r;
                t3.i = cc[ido * (1 + cdim * k)].i + cc[ido * (3 + cdim * k)].i;
                t4.r = cc[ido * (1 + cdim * k)].r - cc[ido * (3 + cdim * k)].r;
                t4.i = cc[ido * (1 + cdim * k)].i - cc[ido * (3 + cdim * k)].i;
            }
            {
                let tmp = t4.r;
                t4.r = -t4.i;
                t4.i = tmp;
            }
            {
                ch[ido * k].r = t2.r + t3.r;
                ch[ido * k].i = t2.i + t3.i;
                ch[ido * (k + l1 * 2)].r = t2.r - t3.r;
                ch[ido * (k + l1 * 2)].i = t2.i - t3.i;
            }
            {
                ch[ido * (k + l1)].r = t1.r + t4.r;
                ch[ido * (k + l1)].i = t1.i + t4.i;
                ch[ido * (k + l1 * 3)].r = t1.r - t4.r;
                ch[ido * (k + l1 * 3)].i = t1.i - t4.i;
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t2.r = cc[ido * cdim * k].r + cc[ido * (2 + cdim * k)].r;
                    t2.i = cc[ido * cdim * k].i + cc[ido * (2 + cdim * k)].i;
                    t1.r = cc[ido * cdim * k].r - cc[ido * (2 + cdim * k)].r;
                    t1.i = cc[ido * cdim * k].i - cc[ido * (2 + cdim * k)].i;
                }
                {
                    t3.r = cc[ido * (1 + cdim * k)].r + cc[ido * (3 + cdim * k)].r;
                    t3.i = cc[ido * (1 + cdim * k)].i + cc[ido * (3 + cdim * k)].i;
                    t4.r = cc[ido * (1 + cdim * k)].r - cc[ido * (3 + cdim * k)].r;
                    t4.i = cc[ido * (1 + cdim * k)].i - cc[ido * (3 + cdim * k)].i;
                }
                {
                    let tmp = t4.r;
                    t4.r = -t4.i;
                    t4.i = tmp;
                }
                {
                    ch[ido * k].r = t2.r + t3.r;
                    ch[ido * k].i = t2.i + t3.i;
                    ch[ido * (k + l1 * 2)].r = t2.r - t3.r;
                    ch[ido * (k + l1 * 2)].i = t2.i - t3.i;
                }
                {
                    ch[ido * (k + l1)].r = t1.r + t4.r;
                    ch[ido * (k + l1)].i = t1.i + t4.i;
                    ch[ido * (k + l1 * 3)].r = t1.r - t4.r;
                    ch[ido * (k + l1 * 3)].i = t1.i - t4.i;
                }
            }
            for i in 1..ido {
                let mut c2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut c3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut c4: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                let cc0: cmplx = cc[i + ido * cdim * k];
                let cc1: cmplx = cc[i + ido * (1 + cdim * k)];
                let cc2: cmplx = cc[i + ido * (2 + cdim * k)];
                let cc3: cmplx = cc[i + ido * (3 + cdim * k)];
                {
                    t2.r = cc0.r + cc2.r;
                    t2.i = cc0.i + cc2.i;
                    t1.r = cc0.r - cc2.r;
                    t1.i = cc0.i - cc2.i;
                }
                {
                    t3.r = cc1.r + cc3.r;
                    t3.i = cc1.i + cc3.i;
                    t4.r = cc1.r - cc3.r;
                    t4.i = cc1.i - cc3.i;
                }
                {
                    let tmp = t4.r;
                    t4.r = -t4.i;
                    t4.i = tmp;
                }
                let wa0: cmplx = wa[i - 1 + 0 * (ido - 1)];
                let wa1: cmplx = wa[i - 1 + ido - 1];
                let wa2: cmplx = wa[i - 1 + 2 * (ido - 1)];
                {
                    ch[i + ido * k].r = t2.r + t3.r;
                    ch[i + ido * k].i = t2.i + t3.i;
                    c3.r = t2.r - t3.r;
                    c3.i = t2.i - t3.i;
                }
                {
                    c2.r = t1.r + t4.r;
                    c2.i = t1.i + t4.i;
                    c4.r = t1.r - t4.r;
                    c4.i = t1.i - t4.i;
                }
                {
                    ch[i + ido * (k + l1)].r = wa0.r * c2.r - wa0.i * c2.i;
                    ch[i + ido * (k + l1)].i = wa0.r * c2.i + wa0.i * c2.r;
                }
                {
                    ch[i + ido * (k + l1 * 2)].r = wa1.r * c3.r - wa1.i * c3.i;
                    ch[i + ido * (k + l1 * 2)].i = wa1.r * c3.i + wa1.i * c3.r;
                }
                {
                    ch[i + ido * (k + l1 * 3)].r = wa2.r * c4.r - wa2.i * c4.i;
                    ch[i + ido * (k + l1 * 3)].i = wa2.r * c4.i + wa2.i * c4.r;
                }
            }
        }
    }
}

fn pass4f(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx]) {
    let cdim: usize = 4;

    if ido == 1 {
        for k in 0..l1 {
            let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
            {
                t2.r = cc[ido * cdim * k].r + cc[ido * (2 + cdim * k)].r;
                t2.i = cc[ido * cdim * k].i + cc[ido * (2 + cdim * k)].i;
                t1.r = cc[ido * cdim * k].r - cc[ido * (2 + cdim * k)].r;
                t1.i = cc[ido * cdim * k].i - cc[ido * (2 + cdim * k)].i;
            }
            {
                t3.r = cc[ido * (1 + cdim * k)].r + cc[ido * (3 + cdim * k)].r;
                t3.i = cc[ido * (1 + cdim * k)].i + cc[ido * (3 + cdim * k)].i;
                t4.r = cc[ido * (1 + cdim * k)].r - cc[ido * (3 + cdim * k)].r;
                t4.i = cc[ido * (1 + cdim * k)].i - cc[ido * (3 + cdim * k)].i;
            }
            {
                let tmp = -t4.r;
                t4.r = t4.i;
                t4.i = tmp;
            }
            {
                ch[ido * k].r = t2.r + t3.r;
                ch[ido * k].i = t2.i + t3.i;
                ch[ido * (k + l1 * 2)].r = t2.r - t3.r;
                ch[ido * (k + l1 * 2)].i = t2.i - t3.i;
            }
            {
                ch[ido * (k + l1)].r = t1.r + t4.r;
                ch[ido * (k + l1)].i = t1.i + t4.i;
                ch[ido * (k + l1 * 3)].r = t1.r - t4.r;
                ch[ido * (k + l1 * 3)].i = t1.i - t4.i;
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t2.r = cc[ido * cdim * k].r + cc[ido * (2 + cdim * k)].r;
                    t2.i = cc[ido * cdim * k].i + cc[ido * (2 + cdim * k)].i;
                    t1.r = cc[ido * cdim * k].r - cc[ido * (2 + cdim * k)].r;
                    t1.i = cc[ido * cdim * k].i - cc[ido * (2 + cdim * k)].i;
                }
                {
                    t3.r = cc[ido * (1 + cdim * k)].r + cc[ido * (3 + cdim * k)].r;
                    t3.i = cc[ido * (1 + cdim * k)].i + cc[ido * (3 + cdim * k)].i;
                    t4.r = cc[ido * (1 + cdim * k)].r - cc[ido * (3 + cdim * k)].r;
                    t4.i = cc[ido * (1 + cdim * k)].i - cc[ido * (3 + cdim * k)].i;
                }
                {
                    let tmp = -t4.r;
                    t4.r = t4.i;
                    t4.i = tmp;
                }
                {
                    ch[ido * k].r = t2.r + t3.r;
                    ch[ido * k].i = t2.i + t3.i;
                    ch[ido * (k + l1 * 2)].r = t2.r - t3.r;
                    ch[ido * (k + l1 * 2)].i = t2.i - t3.i;
                }
                {
                    ch[ido * (k + l1)].r = t1.r + t4.r;
                    ch[ido * (k + l1)].i = t1.i + t4.i;
                    ch[ido * (k + l1 * 3)].r = t1.r - t4.r;
                    ch[ido * (k + l1 * 3)].i = t1.i - t4.i;
                }
            }
            for i in 1..ido {
                let mut c2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut c3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut c4: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                let cc0: cmplx = cc[i + ido * cdim * k];
                let cc1: cmplx = cc[i + ido * (1 + cdim * k)];
                let cc2: cmplx = cc[i + ido * (2 + cdim * k)];
                let cc3: cmplx = cc[i + ido * (3 + cdim * k)];
                {
                    t2.r = cc0.r + cc2.r;
                    t2.i = cc0.i + cc2.i;
                    t1.r = cc0.r - cc2.r;
                    t1.i = cc0.i - cc2.i;
                }
                {
                    t3.r = cc1.r + cc3.r;
                    t3.i = cc1.i + cc3.i;
                    t4.r = cc1.r - cc3.r;
                    t4.i = cc1.i - cc3.i;
                }
                {
                    let tmp = -t4.r;
                    t4.r = t4.i;
                    t4.i = tmp;
                }
                let wa0: cmplx = wa[i - 1 + 0 * (ido - 1)];
                let wa1: cmplx = wa[i - 1 + (ido - 1)];
                let wa2: cmplx = wa[i - 1 + 2 * (ido - 1)];
                {
                    ch[i + ido * k].r = t2.r + t3.r;
                    ch[i + ido * k].i = t2.i + t3.i;
                    c3.r = t2.r - t3.r;
                    c3.i = t2.i - t3.i;
                }
                {
                    c2.r = t1.r + t4.r;
                    c2.i = t1.i + t4.i;
                    c4.r = t1.r - t4.r;
                    c4.i = t1.i - t4.i;
                }
                {
                    ch[i + ido * (k + l1)].r = wa0.r * c2.r + wa0.i * c2.i;
                    ch[i + ido * (k + l1)].i = wa0.r * c2.i - wa0.i * c2.r;
                }
                {
                    ch[i + ido * (k + l1 * 2)].r = wa1.r * c3.r + wa1.i * c3.i;
                    ch[i + ido * (k + l1 * 2)].i = wa1.r * c3.i - wa1.i * c3.r;
                }
                {
                    ch[i + ido * (k + l1 * 3)].r = wa2.r * c4.r + wa2.i * c4.i;
                    ch[i + ido * (k + l1 * 3)].i = wa2.r * c4.i - wa2.i * c4.r;
                }
            }
        }
    }
}

fn pass5b(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx]) {
    let cdim: usize = 5;
    let tw1r: f64 = 0.309_016_994_374_947_45;
    let tw1i: f64 = 0.951_056_516_295_153_5;
    let tw2r: f64 = -0.809_016_994_374_947_5;
    let tw2i: f64 = 0.587_785_252_292_473_1;

    if ido == 1 {
        for k in 0..l1 {
            let t0: cmplx = cc[ido * cdim * k];
            let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
            {
                t1.r = cc[ido * (1 + cdim * k)].r + cc[ido * (4 + cdim * k)].r;
                t1.i = cc[ido * (1 + cdim * k)].i + cc[ido * (4 + cdim * k)].i;
                t4.r = cc[ido * (1 + cdim * k)].r - cc[ido * (4 + cdim * k)].r;
                t4.i = cc[ido * (1 + cdim * k)].i - cc[ido * (4 + cdim * k)].i;
            }
            {
                t2.r = cc[ido * (2 + cdim * k)].r + cc[ido * (3 + cdim * k)].r;
                t2.i = cc[ido * (2 + cdim * k)].i + cc[ido * (3 + cdim * k)].i;
                t3.r = cc[ido * (2 + cdim * k)].r - cc[ido * (3 + cdim * k)].r;
                t3.i = cc[ido * (2 + cdim * k)].i - cc[ido * (3 + cdim * k)].i;
            }
            ch[ido * k].r = t0.r + t1.r + t2.r;
            ch[ido * k].i = t0.i + t1.i + t2.i;
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t0.r + tw1r * t1.r + tw2r * t2.r;
                ca.i = t0.i + tw1r * t1.i + tw2r * t2.i;
                cb.i = tw1i * t4.r + tw2i * t3.r;
                cb.r = -(tw1i * t4.i + tw2i * t3.i);
                {
                    ch[ido * (k + l1)].r = ca.r + cb.r;
                    ch[ido * (k + l1)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 4)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 4)].i = ca.i - cb.i;
                }
            }
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t0.r + tw2r * t1.r + tw1r * t2.r;
                ca.i = t0.i + tw2r * t1.i + tw1r * t2.i;
                cb.i = tw2i * t4.r - tw1i * t3.r;
                cb.r = -(tw2i * t4.i - tw1i * t3.i);
                {
                    ch[ido * (k + l1 * 2)].r = ca.r + cb.r;
                    ch[ido * (k + l1 * 2)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 3)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 3)].i = ca.i - cb.i;
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let t0: cmplx = cc[ido * cdim * k];
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t1.r = cc[ido * (1 + cdim * k)].r + cc[ido * (4 + cdim * k)].r;
                    t1.i = cc[ido * (1 + cdim * k)].i + cc[ido * (4 + cdim * k)].i;
                    t4.r = cc[ido * (1 + cdim * k)].r - cc[ido * (4 + cdim * k)].r;
                    t4.i = cc[ido * (1 + cdim * k)].i - cc[ido * (4 + cdim * k)].i;
                }
                {
                    t2.r = cc[ido * (2 + cdim * k)].r + cc[ido * (3 + cdim * k)].r;
                    t2.i = cc[ido * (2 + cdim * k)].i + cc[ido * (3 + cdim * k)].i;
                    t3.r = cc[ido * (2 + cdim * k)].r - cc[ido * (3 + cdim * k)].r;
                    t3.i = cc[ido * (2 + cdim * k)].i - cc[ido * (3 + cdim * k)].i;
                }
                ch[ido * k].r = t0.r + t1.r + t2.r;
                ch[ido * k].i = t0.i + t1.i + t2.i;
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw1r * t1.r + tw2r * t2.r;
                    ca.i = t0.i + tw1r * t1.i + tw2r * t2.i;
                    cb.i = tw1i * t4.r + tw2i * t3.r;
                    cb.r = -(tw1i * t4.i + tw2i * t3.i);
                    {
                        ch[ido * (k + l1)].r = ca.r + cb.r;
                        ch[ido * (k + l1)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 4)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 4)].i = ca.i - cb.i;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw2r * t1.r + tw1r * t2.r;
                    ca.i = t0.i + tw2r * t1.i + tw1r * t2.i;
                    cb.i = tw2i * t4.r - tw1i * t3.r;
                    cb.r = -(tw2i * t4.i - tw1i * t3.i);
                    {
                        ch[ido * (k + l1 * 2)].r = ca.r + cb.r;
                        ch[ido * (k + l1 * 2)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 3)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 3)].i = ca.i - cb.i;
                    }
                }
            }
            for i in 1..ido {
                let t0: cmplx = cc[i + ido * cdim * k];
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t1.r = cc[i + ido * (1 + cdim * k)].r + cc[i + ido * (4 + cdim * k)].r;
                    t1.i = cc[i + ido * (1 + cdim * k)].i + cc[i + ido * (4 + cdim * k)].i;
                    t4.r = cc[i + ido * (1 + cdim * k)].r - cc[i + ido * (4 + cdim * k)].r;
                    t4.i = cc[i + ido * (1 + cdim * k)].i - cc[i + ido * (4 + cdim * k)].i;
                }
                {
                    t2.r = cc[i + ido * (2 + cdim * k)].r + cc[i + ido * (3 + cdim * k)].r;
                    t2.i = cc[i + ido * (2 + cdim * k)].i + cc[i + ido * (3 + cdim * k)].i;
                    t3.r = cc[i + ido * (2 + cdim * k)].r - cc[i + ido * (3 + cdim * k)].r;
                    t3.i = cc[i + ido * (2 + cdim * k)].i - cc[i + ido * (3 + cdim * k)].i;
                }
                ch[i + ido * k].r = t0.r + t1.r + t2.r;
                ch[i + ido * k].i = t0.i + t1.i + t2.i;
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw1r * t1.r + tw2r * t2.r;
                    ca.i = t0.i + tw1r * t1.i + tw2r * t2.i;
                    cb.i = tw1i * t4.r + tw2i * t3.r;
                    cb.r = -(tw1i * t4.i + tw2i * t3.i);
                    {
                        da.r = ca.r + cb.r;
                        da.i = ca.i + cb.i;
                        db.r = ca.r - cb.r;
                        db.i = ca.i - cb.i;
                    }
                    {
                        ch[i + ido * (k + l1)].r = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.r
                            - wa[i - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1)].i = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.i
                            + wa[i - 1 + (1 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 4)].r = wa[i - 1 + (4 - 1) * (ido - 1)].r * db.r
                            - wa[i - 1 + (4 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 4)].i = wa[i - 1 + (4 - 1) * (ido - 1)].r * db.i
                            + wa[i - 1 + (4 - 1) * (ido - 1)].i * db.r;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw2r * t1.r + tw1r * t2.r;
                    ca.i = t0.i + tw2r * t1.i + tw1r * t2.i;
                    cb.i = tw2i * t4.r - tw1i * t3.r;
                    cb.r = -(tw2i * t4.i - tw1i * t3.i);
                    {
                        da.r = ca.r + cb.r;
                        da.i = ca.i + cb.i;
                        db.r = ca.r - cb.r;
                        db.i = ca.i - cb.i;
                    }
                    {
                        ch[i + ido * (k + l1 * 2)].r =
                            wa[i - 1 + (ido - 1)].r * da.r - wa[i - 1 + (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1 * 2)].i =
                            wa[i - 1 + (ido - 1)].r * da.i + wa[i - 1 + (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 3)].r = wa[i - 1 + (3 - 1) * (ido - 1)].r * db.r
                            - wa[i - 1 + (3 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 3)].i = wa[i - 1 + (3 - 1) * (ido - 1)].r * db.i
                            + wa[i - 1 + (3 - 1) * (ido - 1)].i * db.r;
                    }
                }
            }
        }
    }
}

fn pass5f(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx]) {
    let cdim: usize = 5;
    let tw1r: f64 = 0.309_016_994_374_947_45;
    let tw1i: f64 = -0.951_056_516_295_153_5;
    let tw2r: f64 = -0.809_016_994_374_947_5;
    let tw2i: f64 = -0.587_785_252_292_473_1;

    if ido == 1 {
        for k in 0..l1 {
            let t0: cmplx = cc[ido * cdim * k];
            let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
            {
                t1.r = cc[ido * (1 + cdim * k)].r + cc[ido * (4 + cdim * k)].r;
                t1.i = cc[ido * (1 + cdim * k)].i + cc[ido * (4 + cdim * k)].i;
                t4.r = cc[ido * (1 + cdim * k)].r - cc[ido * (4 + cdim * k)].r;
                t4.i = cc[ido * (1 + cdim * k)].i - cc[ido * (4 + cdim * k)].i;
            }
            {
                t2.r = cc[ido * (2 + cdim * k)].r + cc[ido * (3 + cdim * k)].r;
                t2.i = cc[ido * (2 + cdim * k)].i + cc[ido * (3 + cdim * k)].i;
                t3.r = cc[ido * (2 + cdim * k)].r - cc[ido * (3 + cdim * k)].r;
                t3.i = cc[ido * (2 + cdim * k)].i - cc[ido * (3 + cdim * k)].i;
            }
            ch[ido * k].r = t0.r + t1.r + t2.r;
            ch[ido * k].i = t0.i + t1.i + t2.i;
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t0.r + tw1r * t1.r + tw2r * t2.r;
                ca.i = t0.i + tw1r * t1.i + tw2r * t2.i;
                cb.i = tw1i * t4.r + tw2i * t3.r;
                cb.r = -(tw1i * t4.i + tw2i * t3.i);
                {
                    ch[ido * (k + l1)].r = ca.r + cb.r;
                    ch[ido * (k + l1)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 4)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 4)].i = ca.i - cb.i;
                }
            }
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };

                ca.r = t0.r + tw2r * t1.r + tw1r * t2.r;
                ca.i = t0.i + tw2r * t1.i + tw1r * t2.i;
                cb.i = tw2i * t4.r - tw1i * t3.r;
                cb.r = -(tw2i * t4.i - tw1i * t3.i);
                {
                    ch[ido * (k + l1 * 2)].r = ca.r + cb.r;
                    ch[ido * (k + l1 * 2)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 3)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 3)].i = ca.i - cb.i;
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let t0: cmplx = cc[ido * cdim * k];
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t1.r = cc[ido * (1 + cdim * k)].r + cc[ido * (4 + cdim * k)].r;
                    t1.i = cc[ido * (1 + cdim * k)].i + cc[ido * (4 + cdim * k)].i;
                    t4.r = cc[ido * (1 + cdim * k)].r - cc[ido * (4 + cdim * k)].r;
                    t4.i = cc[ido * (1 + cdim * k)].i - cc[ido * (4 + cdim * k)].i;
                }
                {
                    t2.r = cc[ido * (2 + cdim * k)].r + cc[ido * (3 + cdim * k)].r;
                    t2.i = cc[ido * (2 + cdim * k)].i + cc[ido * (3 + cdim * k)].i;
                    t3.r = cc[ido * (2 + cdim * k)].r - cc[ido * (3 + cdim * k)].r;
                    t3.i = cc[ido * (2 + cdim * k)].i - cc[ido * (3 + cdim * k)].i;
                }
                ch[ido * k].r = t0.r + t1.r + t2.r;
                ch[ido * k].i = t0.i + t1.i + t2.i;
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw1r * t1.r + tw2r * t2.r;
                    ca.i = t0.i + tw1r * t1.i + tw2r * t2.i;
                    cb.i = tw1i * t4.r + tw2i * t3.r;
                    cb.r = -(tw1i * t4.i + tw2i * t3.i);
                    {
                        ch[ido * (k + l1)].r = ca.r + cb.r;
                        ch[ido * (k + l1)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 4)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 4)].i = ca.i - cb.i;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw2r * t1.r + tw1r * t2.r;
                    ca.i = t0.i + tw2r * t1.i + tw1r * t2.i;
                    cb.i = tw2i * t4.r - tw1i * t3.r;
                    cb.r = -(tw2i * t4.i - tw1i * t3.i);
                    {
                        ch[ido * (k + l1 * 2)].r = ca.r + cb.r;
                        ch[ido * (k + l1 * 2)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 3)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 3)].i = ca.i - cb.i;
                    }
                }
            }
            for i in 1..ido {
                let t0: cmplx = cc[(i) + ido * cdim * (k)];
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t1.r = cc[i + ido * (1 + cdim * k)].r + cc[i + ido * (4 + cdim * k)].r;
                    t1.i = cc[i + ido * (1 + cdim * k)].i + cc[i + ido * (4 + cdim * k)].i;
                    t4.r = cc[i + ido * (1 + cdim * k)].r - cc[i + ido * (4 + cdim * k)].r;
                    t4.i = cc[i + ido * (1 + cdim * k)].i - cc[i + ido * (4 + cdim * k)].i;
                }
                {
                    t2.r = cc[i + ido * (2 + cdim * k)].r + cc[i + ido * (3 + cdim * k)].r;
                    t2.i = cc[i + ido * (2 + cdim * k)].i + cc[i + ido * (3 + cdim * k)].i;
                    t3.r = cc[i + ido * (2 + cdim * k)].r - cc[i + ido * (3 + cdim * k)].r;
                    t3.i = cc[i + ido * (2 + cdim * k)].i - cc[i + ido * (3 + cdim * k)].i;
                }
                ch[i + ido * k].r = t0.r + t1.r + t2.r;
                ch[i + ido * k].i = t0.i + t1.i + t2.i;
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw1r * t1.r + tw2r * t2.r;
                    ca.i = t0.i + tw1r * t1.i + tw2r * t2.i;
                    cb.i = tw1i * t4.r + tw2i * t3.r;
                    cb.r = -(tw1i * t4.i + tw2i * t3.i);
                    {
                        da.r = ca.r + cb.r;
                        da.i = ca.i + cb.i;
                        db.r = ca.r - cb.r;
                        db.i = ca.i - cb.i;
                    }
                    {
                        ch[i + ido * (k + l1)].r = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.r
                            + wa[i - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1)].i = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.i
                            - wa[i - 1 + (1 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 4)].r = wa[i - 1 + (4 - 1) * (ido - 1)].r * db.r
                            + wa[i - 1 + (4 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 4)].i = wa[i - 1 + (4 - 1) * (ido - 1)].r * db.i
                            - wa[i - 1 + (4 - 1) * (ido - 1)].i * db.r;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t0.r + tw2r * t1.r + tw1r * t2.r;
                    ca.i = t0.i + tw2r * t1.i + tw1r * t2.i;
                    cb.i = tw2i * t4.r - tw1i * t3.r;
                    cb.r = -(tw2i * t4.i - tw1i * t3.i);
                    {
                        da.r = ca.r + cb.r;
                        da.i = ca.i + cb.i;
                        db.r = ca.r - cb.r;
                        db.i = ca.i - cb.i;
                    }
                    {
                        ch[i + ido * (k + l1 * 2)].r =
                            wa[(i) - 1 + (ido - 1)].r * da.r + wa[i - 1 + (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1 * 2)].i =
                            wa[i - 1 + (ido - 1)].r * da.i - wa[i - 1 + (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 3)].r = wa[i - 1 + (3 - 1) * (ido - 1)].r * db.r
                            + wa[i - 1 + (3 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 3)].i = wa[i - 1 + (3 - 1) * (ido - 1)].r * db.i
                            - wa[i - 1 + (3 - 1) * (ido - 1)].i * db.r;
                    }
                }
            }
        }
    }
}

fn pass7(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx], sign: i64) {
    let cdim: usize = 7;
    let tw1r: f64 = 0.623_489_801_858_733_5;
    let tw1i: f64 = (sign as f64) * 0.781_831_482_468_029_8;
    let tw2r = -0.222_520_933_956_314_4;
    let tw2i: f64 = (sign as f64) * 0.974_927_912_181_823_6;
    let tw3r: f64 = -0.900_968_867_902_419_1;
    let tw3i: f64 = (sign as f64) * 0.433_883_739_117_558_1;

    if ido == 1 {
        for k in 0..l1 {
            let t1: cmplx = cc[ido * cdim * k];
            let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t5: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t6: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t7: cmplx = cmplx { r: 0.0, i: 0.0 };
            {
                t2.r = cc[ido * (1 + cdim * k)].r + cc[ido * (6 + cdim * k)].r;
                t2.i = cc[ido * (1 + cdim * k)].i + cc[ido * (6 + cdim * k)].i;
                t7.r = cc[ido * (1 + cdim * k)].r - cc[ido * (6 + cdim * k)].r;
                t7.i = cc[ido * (1 + cdim * k)].i - cc[ido * (6 + cdim * k)].i;
            }
            {
                t3.r = cc[ido * (2 + cdim * k)].r + cc[ido * (5 + cdim * k)].r;
                t3.i = cc[ido * (2 + cdim * k)].i + cc[ido * (5 + cdim * k)].i;
                t6.r = cc[ido * (2 + cdim * k)].r - cc[ido * (5 + cdim * k)].r;
                t6.i = cc[ido * (2 + cdim * k)].i - cc[ido * (5 + cdim * k)].i;
            }
            {
                t4.r = cc[ido * (3 + cdim * k)].r + cc[ido * (4 + cdim * k)].r;
                t4.i = cc[ido * (3 + cdim * k)].i + cc[ido * (4 + cdim * k)].i;
                t5.r = cc[ido * (3 + cdim * k)].r - cc[ido * (4 + cdim * k)].r;
                t5.i = cc[ido * (3 + cdim * k)].i - cc[ido * (4 + cdim * k)].i;
            }
            ch[ido * k].r = t1.r + t2.r + t3.r + t4.r;
            ch[ido * k].i = t1.i + t2.i + t3.i + t4.i;
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r;
                ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i;
                cb.i = tw1i * t7.r + tw2i * t6.r + tw3i * t5.r;
                cb.r = -(tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                {
                    ch[ido * (k + l1)].r = ca.r + cb.r;
                    ch[ido * (k + l1)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 6)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 6)].i = ca.i - cb.i;
                }
            }
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r;
                ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i;
                cb.i = tw2i * t7.r - tw3i * t6.r - tw1i * t5.r;
                cb.r = -(tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                {
                    ch[ido * (k + l1 * 2)].r = ca.r + cb.r;
                    ch[ido * (k + l1 * 2)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 5)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 5)].i = ca.i - cb.i;
                }
            }
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r;
                ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i;
                cb.i = tw3i * t7.r - tw1i * t6.r + tw2i * t5.r;
                cb.r = -(tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                {
                    ch[ido * (k + l1 * 3)].r = ca.r + cb.r;
                    ch[ido * (k + l1 * 3)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 4)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 4)].i = ca.i - cb.i;
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let t1: cmplx = cc[ido * cdim * k];
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t5: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t6: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t7: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t2.r = cc[ido * (1 + cdim * k)].r + cc[ido * (6 + cdim * k)].r;
                    t2.i = cc[ido * (1 + cdim * k)].i + cc[ido * (6 + cdim * k)].i;
                    t7.r = cc[ido * (1 + cdim * k)].r - cc[ido * (6 + cdim * k)].r;
                    t7.i = cc[ido * (1 + cdim * k)].i - cc[ido * (6 + cdim * k)].i;
                }
                {
                    t3.r = cc[ido * (2 + cdim * k)].r + cc[ido * (5 + cdim * k)].r;
                    t3.i = cc[ido * (2 + cdim * k)].i + cc[ido * (5 + cdim * k)].i;
                    t6.r = cc[ido * (2 + cdim * k)].r - cc[ido * (5 + cdim * k)].r;
                    t6.i = cc[ido * (2 + cdim * k)].i - cc[ido * (5 + cdim * k)].i;
                }
                {
                    t4.r = cc[ido * (3 + cdim * k)].r + cc[ido * (4 + cdim * k)].r;
                    t4.i = cc[ido * (3 + cdim * k)].i + cc[ido * (4 + cdim * k)].i;
                    t5.r = cc[ido * (3 + cdim * k)].r - cc[ido * (4 + cdim * k)].r;
                    t5.i = cc[ido * (3 + cdim * k)].i - cc[ido * (4 + cdim * k)].i;
                }
                ch[ido * k].r = t1.r + t2.r + t3.r + t4.r;
                ch[ido * k].i = t1.i + t2.i + t3.i + t4.i;

                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r;
                    ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i;
                    cb.i = tw1i * t7.r + tw2i * t6.r + tw3i * t5.r;
                    cb.r = -(tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                    {
                        ch[ido * (k + l1)].r = ca.r + cb.r;
                        ch[ido * (k + l1)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 6)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 6)].i = ca.i - cb.i;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r;
                    ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i;
                    cb.i = tw2i * t7.r - tw3i * t6.r - tw1i * t5.r;
                    cb.r = -(tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                    {
                        ch[ido * (k + l1 * 2)].r = ca.r + cb.r;
                        ch[ido * (k + l1 * 2)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 5)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 5)].i = ca.i - cb.i;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r;
                    ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i;
                    cb.i = tw3i * t7.r - tw1i * t6.r + tw2i * t5.r;
                    cb.r = -(tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                    {
                        ch[ido * (k + l1 * 3)].r = ca.r + cb.r;
                        ch[ido * (k + l1 * 3)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 4)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 4)].i = ca.i - cb.i;
                    }
                }
            }
            for i in 1..ido {
                let t1: cmplx = cc[i + ido * cdim * k];
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t5: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t6: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t7: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t2.r = cc[i + ido * (1 + cdim * k)].r + cc[i + ido * (6 + cdim * k)].r;
                    t2.i = cc[i + ido * (1 + cdim * k)].i + cc[i + ido * (6 + cdim * k)].i;
                    t7.r = cc[i + ido * (1 + cdim * k)].r - cc[i + ido * (6 + cdim * k)].r;
                    t7.i = cc[i + ido * (1 + cdim * k)].i - cc[i + ido * (6 + cdim * k)].i;
                }
                {
                    t3.r = cc[i + ido * (2 + cdim * k)].r + cc[i + ido * (5 + cdim * k)].r;
                    t3.i = cc[i + ido * (2 + cdim * k)].i + cc[i + ido * (5 + cdim * k)].i;
                    t6.r = cc[i + ido * (2 + cdim * k)].r - cc[i + ido * (5 + cdim * k)].r;
                    t6.i = cc[i + ido * (2 + cdim * k)].i - cc[i + ido * (5 + cdim * k)].i;
                }
                {
                    t4.r = cc[i + ido * (3 + cdim * k)].r + cc[i + ido * (4 + cdim * k)].r;
                    t4.i = cc[i + ido * (3 + cdim * k)].i + cc[i + ido * (4 + cdim * k)].i;
                    t5.r = cc[i + ido * (3 + cdim * k)].r - cc[i + ido * (4 + cdim * k)].r;
                    t5.i = cc[i + ido * (3 + cdim * k)].i - cc[i + ido * (4 + cdim * k)].i;
                }
                ch[i + ido * k].r = t1.r + t2.r + t3.r + t4.r;
                ch[i + ido * k].i = t1.i + t2.i + t3.i + t4.i;
                {
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                        let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                        ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r;
                        ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i;
                        cb.i = tw1i * t7.r + tw2i * t6.r + tw3i * t5.r;
                        cb.r = -(tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                        {
                            da.r = ca.r + cb.r;
                            da.i = ca.i + cb.i;
                            db.r = ca.r - cb.r;
                            db.i = ca.i - cb.i;
                        }
                    }
                    {
                        ch[i + ido * (k + l1)].r = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.r
                            - (sign as f64) * wa[i - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1)].i = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.i
                            + (sign as f64) * wa[i - 1 + (1 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 6)].r = wa[i - 1 + (6 - 1) * (ido - 1)].r * db.r
                            - (sign as f64) * wa[i - 1 + (6 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 6)].i = wa[i - 1 + (6 - 1) * (ido - 1)].r * db.i
                            + (sign as f64) * wa[i - 1 + (6 - 1) * (ido - 1)].i * db.r;
                    }
                }
                {
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                        let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                        ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r;
                        ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i;
                        cb.i = tw2i * t7.r - tw3i * t6.r - tw1i * t5.r;
                        cb.r = -(tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                        {
                            da.r = ca.r + cb.r;
                            da.i = ca.i + cb.i;
                            db.r = ca.r - cb.r;
                            db.i = ca.i - cb.i;
                        }
                    }
                    {
                        ch[i + ido * (k + l1 * 2)].r = wa[i - 1 + (ido - 1)].r * da.r
                            - (sign as f64) * wa[i - 1 + (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1 * 2)].i = wa[i - 1 + (ido - 1)].r * da.i
                            + (sign as f64) * wa[i - 1 + (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 5)].r = wa[i - 1 + (5 - 1) * (ido - 1)].r * db.r
                            - (sign as f64) * wa[i - 1 + (5 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 5)].i = wa[i - 1 + (5 - 1) * (ido - 1)].r * db.i
                            + (sign as f64) * wa[i - 1 + (5 - 1) * (ido - 1)].i * db.r;
                    }
                }
                {
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                        let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                        ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r;
                        ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i;
                        cb.i = tw3i * t7.r - tw1i * t6.r + tw2i * t5.r;
                        cb.r = -(tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                        {
                            da.r = ca.r + cb.r;
                            da.i = ca.i + cb.i;
                            db.r = ca.r - cb.r;
                            db.i = ca.i - cb.i;
                        }
                    }
                    {
                        ch[i + ido * (k + l1 * 3)].r = wa[i - 1 + (3 - 1) * (ido - 1)].r * da.r
                            - (sign as f64) * wa[i - 1 + (3 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1 * 3)].i = wa[i - 1 + (3 - 1) * (ido - 1)].r * da.i
                            + (sign as f64) * wa[i - 1 + (3 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 4)].r = wa[i - 1 + (4 - 1) * (ido - 1)].r * db.r
                            - (sign as f64) * wa[i - 1 + (4 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 4)].i = wa[i - 1 + (4 - 1) * (ido - 1)].r * db.i
                            + (sign as f64) * wa[i - 1 + (4 - 1) * (ido - 1)].i * db.r;
                    }
                }
            }
        }
    }
}

fn pass11(ido: usize, l1: usize, cc: &[cmplx], ch: &mut [cmplx], wa: &[cmplx], sign: i64) {
    let cdim: usize = 11;
    let tw1r: f64 = 0.841_253_532_831_181_2;
    let tw1i: f64 = (sign as f64) * 0.540_640_817_455_597_6;
    let tw2r: f64 = 0.415_415_013_001_886_44;
    let tw2i: f64 = (sign as f64) * 0.909_631_995_354_518_3;
    let tw3r: f64 = -0.142_314_838_273_285_14;
    let tw3i: f64 = (sign as f64) * 0.989_821_441_880_932_7;
    let tw4r: f64 = -0.654_860_733_945_285_1;
    let tw4i: f64 = (sign as f64) * 0.755_749_574_354_258_3;
    let tw5r: f64 = -0.959_492_973_614_497_4;
    let tw5i: f64 = (sign as f64) * 0.281_732_556_841_429_67;

    if ido == 1 {
        for k in 0..l1 {
            let t1: cmplx = cc[ido * cdim * k];
            let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t5: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t6: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t7: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t8: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t9: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t10: cmplx = cmplx { r: 0.0, i: 0.0 };
            let mut t11: cmplx = cmplx { r: 0.0, i: 0.0 };
            {
                t2.r = &cc[ido * (1 + cdim * k)].r + cc[ido * (10 + cdim * k)].r;
                t2.i = &cc[ido * (1 + cdim * k)].i + cc[ido * (10 + cdim * k)].i;
                t11.r = &cc[ido * (1 + cdim * k)].r - cc[ido * (10 + cdim * k)].r;
                t11.i = &cc[ido * (1 + cdim * k)].i - cc[ido * (10 + cdim * k)].i;
            }
            {
                t3.r = &cc[ido * (2 + cdim * k)].r + cc[ido * (9 + cdim * k)].r;
                t3.i = &cc[ido * (2 + cdim * k)].i + cc[ido * (9 + cdim * k)].i;
                t10.r = &cc[ido * (2 + cdim * k)].r - cc[ido * (9 + cdim * k)].r;
                t10.i = &cc[ido * (2 + cdim * k)].i - cc[ido * (9 + cdim * k)].i;
            }
            {
                t4.r = &cc[ido * (3 + cdim * k)].r + cc[ido * (8 + cdim * k)].r;
                t4.i = &cc[ido * (3 + cdim * k)].i + cc[ido * (8 + cdim * k)].i;
                t9.r = &cc[ido * (3 + cdim * k)].r - cc[ido * (8 + cdim * k)].r;
                t9.i = &cc[ido * (3 + cdim * k)].i - cc[ido * (8 + cdim * k)].i;
            }
            {
                t5.r = &cc[ido * (4 + cdim * k)].r + cc[ido * (7 + cdim * k)].r;
                t5.i = &cc[ido * (4 + cdim * k)].i + cc[ido * (7 + cdim * k)].i;
                t8.r = &cc[ido * (4 + cdim * k)].r - cc[ido * (7 + cdim * k)].r;
                t8.i = &cc[ido * (4 + cdim * k)].i - cc[ido * (7 + cdim * k)].i;
            }
            {
                t6.r = cc[ido * (5 + cdim * k)].r + cc[ido * (6 + cdim * k)].r;
                t6.i = cc[ido * (5 + cdim * k)].i + cc[ido * (6 + cdim * k)].i;
                t7.r = cc[ido * (5 + cdim * k)].r - cc[ido * (6 + cdim * k)].r;
                t7.i = cc[ido * (5 + cdim * k)].i - cc[ido * (6 + cdim * k)].i;
            }
            ch[ido * k].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r;
            ch[ido * k].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r;
                ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i;
                cb.i = tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r;
                cb.r = -(tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i);
                {
                    ch[ido * (k + l1)].r = ca.r + cb.r;
                    ch[ido * (k + l1)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 10)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 10)].i = ca.i - cb.i;
                }
            }
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r;
                ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i;
                cb.i = tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r;
                cb.r = -(tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i);
                {
                    ch[ido * (k + l1 * 2)].r = ca.r + cb.r;
                    ch[ido * (k + l1 * 2)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 9)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 9)].i = ca.i - cb.i;
                }
            }
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r;
                ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i;
                cb.i = tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r;
                cb.r = -(tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i);
                {
                    ch[ido * (k + l1 * 3)].r = ca.r + cb.r;
                    ch[ido * (k + l1 * 3)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 8)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 8)].i = ca.i - cb.i;
                }
            }
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r;
                ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i;
                cb.i = tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r;
                cb.r = -(tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i);
                {
                    ch[ido * (k + l1 * 4)].r = ca.r + cb.r;
                    ch[ido * (k + l1 * 4)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 7)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 7)].i = ca.i - cb.i;
                }
            }
            {
                let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r;
                ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i;
                cb.i = tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r;
                cb.r = -(tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i);
                {
                    ch[ido * (k + l1 * 5)].r = ca.r + cb.r;
                    ch[ido * (k + l1 * 5)].i = ca.i + cb.i;
                    ch[ido * (k + l1 * 6)].r = ca.r - cb.r;
                    ch[ido * (k + l1 * 6)].i = ca.i - cb.i;
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t1: cmplx = cmplx { r: 0.0, i: 0.0 };
                t1.r = cc[ido * cdim * k].r;
                t1.i = cc[ido * cdim * k].i;
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t5: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t6: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t7: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t8: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t9: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t10: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t11: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t2.r = cc[ido * (1 + cdim * k)].r + cc[ido * (10 + cdim * k)].r;
                    t2.i = cc[ido * (1 + cdim * k)].i + cc[ido * (10 + cdim * k)].i;
                    t11.r = cc[ido * (1 + cdim * k)].r - cc[ido * (10 + cdim * k)].r;
                    t11.i = cc[ido * (1 + cdim * k)].i - cc[ido * (10 + cdim * k)].i;
                }
                {
                    t3.r = cc[ido * (2 + cdim * k)].r + cc[ido * (9 + cdim * k)].r;
                    t3.i = cc[ido * (2 + cdim * k)].i + cc[ido * (9 + cdim * k)].i;
                    t10.r = cc[ido * (2 + cdim * k)].r - cc[ido * (9 + cdim * k)].r;
                    t10.i = cc[ido * (2 + cdim * k)].i - cc[ido * (9 + cdim * k)].i;
                }
                {
                    t4.r = cc[ido * (3 + cdim * k)].r + cc[ido * (8 + cdim * k)].r;
                    t4.i = cc[ido * (3 + cdim * k)].i + cc[ido * (8 + cdim * k)].i;
                    t9.r = cc[ido * (3 + cdim * k)].r - cc[ido * (8 + cdim * k)].r;
                    t9.i = cc[ido * (3 + cdim * k)].i - cc[ido * (8 + cdim * k)].i;
                }
                {
                    t5.r = cc[ido * (4 + cdim * k)].r + cc[ido * (7 + cdim * k)].r;
                    t5.i = cc[ido * (4 + cdim * k)].i + cc[ido * (7 + cdim * k)].i;
                    t8.r = cc[ido * (4 + cdim * k)].r - cc[ido * (7 + cdim * k)].r;
                    t8.i = cc[ido * (4 + cdim * k)].i - cc[ido * (7 + cdim * k)].i;
                }
                {
                    t6.r = cc[ido * (5 + cdim * k)].r + cc[ido * (6 + cdim * k)].r;
                    t6.i = cc[ido * (5 + cdim * k)].i + cc[ido * (6 + cdim * k)].i;
                    t7.r = cc[ido * (5 + cdim * k)].r - cc[ido * (6 + cdim * k)].r;
                    t7.i = cc[ido * (5 + cdim * k)].i - cc[ido * (6 + cdim * k)].i;
                }
                ch[ido * k].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r;
                ch[ido * k].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r =
                        t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r;
                    ca.i =
                        t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i;
                    cb.i = tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r;
                    cb.r = -(tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i);
                    {
                        ch[ido * (k + l1)].r = ca.r + cb.r;
                        ch[ido * (k + l1)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 10)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 10)].i = ca.i - cb.i;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r =
                        t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r;
                    ca.i =
                        t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i;
                    cb.i = tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r;
                    cb.r = -(tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i);
                    {
                        ch[ido * (k + l1 * 2)].r = ca.r + cb.r;
                        ch[ido * (k + l1 * 2)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 9)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 9)].i = ca.i - cb.i;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r =
                        t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r;
                    ca.i =
                        t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i;
                    cb.i = tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r;
                    cb.r = -(tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i);
                    {
                        ch[ido * (k + l1 * 3)].r = ca.r + cb.r;
                        ch[ido * (k + l1 * 3)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 8)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 8)].i = ca.i - cb.i;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r =
                        t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r;
                    ca.i =
                        t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i;
                    cb.i = tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r;
                    cb.r = -(tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i);
                    {
                        ch[ido * (k + l1 * 4)].r = ca.r + cb.r;
                        ch[ido * (k + l1 * 4)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 7)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 7)].i = ca.i - cb.i;
                    }
                }
                {
                    let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                    ca.r =
                        t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r;
                    ca.i =
                        t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i;
                    cb.i = tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r;
                    cb.r = -(tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i);
                    {
                        ch[ido * (k + l1 * 5)].r = ca.r + cb.r;
                        ch[ido * (k + l1 * 5)].i = ca.i + cb.i;
                        ch[ido * (k + l1 * 6)].r = ca.r - cb.r;
                        ch[ido * (k + l1 * 6)].i = ca.i - cb.i;
                    }
                }
            }
            for i in 1..ido {
                let t1: cmplx = cc[i + ido * cdim * k];
                let mut t2: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t3: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t4: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t5: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t6: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t7: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t8: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t9: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t10: cmplx = cmplx { r: 0.0, i: 0.0 };
                let mut t11: cmplx = cmplx { r: 0.0, i: 0.0 };
                {
                    t2.r = cc[i + ido * (1 + cdim * k)].r + cc[i + ido * (10 + cdim * k)].r;
                    t2.i = cc[i + ido * (1 + cdim * k)].i + cc[i + ido * (10 + cdim * k)].i;
                    t11.r = cc[i + ido * (1 + cdim * k)].r - cc[i + ido * (10 + cdim * k)].r;
                    t11.i = cc[i + ido * (1 + cdim * k)].i - cc[i + ido * (10 + cdim * k)].i;
                }
                {
                    t3.r = cc[i + ido * (2 + cdim * k)].r + cc[i + ido * (9 + cdim * k)].r;
                    t3.i = cc[i + ido * (2 + cdim * k)].i + cc[i + ido * (9 + cdim * k)].i;
                    t10.r = cc[i + ido * (2 + cdim * k)].r - cc[i + ido * (9 + cdim * k)].r;
                    t10.i = cc[i + ido * (2 + cdim * k)].i - cc[i + ido * (9 + cdim * k)].i;
                }
                {
                    t4.r = cc[i + ido * (3 + cdim * k)].r + cc[i + ido * (8 + cdim * k)].r;
                    t4.i = cc[i + ido * (3 + cdim * k)].i + cc[i + ido * (8 + cdim * k)].i;
                    t9.r = cc[i + ido * (3 + cdim * k)].r - cc[i + ido * (8 + cdim * k)].r;
                    t9.i = cc[i + ido * (3 + cdim * k)].i - cc[i + ido * (8 + cdim * k)].i;
                }
                {
                    t5.r = cc[i + ido * (4 + cdim * k)].r + cc[i + ido * (7 + cdim * k)].r;
                    t5.i = cc[i + ido * (4 + cdim * k)].i + cc[i + ido * (7 + cdim * k)].i;
                    t8.r = cc[i + ido * (4 + cdim * k)].r - cc[i + ido * (7 + cdim * k)].r;
                    t8.i = cc[i + ido * (4 + cdim * k)].i - cc[i + ido * (7 + cdim * k)].i;
                }
                {
                    t6.r = cc[i + ido * (5 + cdim * k)].r + cc[i + ido * (6 + cdim * k)].r;
                    t6.i = cc[i + ido * (5 + cdim * k)].i + cc[i + ido * (6 + cdim * k)].i;
                    t7.r = cc[i + ido * (5 + cdim * k)].r - cc[i + ido * (6 + cdim * k)].r;
                    t7.i = cc[i + ido * (5 + cdim * k)].i - cc[i + ido * (6 + cdim * k)].i;
                }
                ch[i + ido * k].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r;
                ch[i + ido * k].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
                {
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                        let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                        ca.r = t1.r
                            + tw1r * t2.r
                            + tw2r * t3.r
                            + tw3r * t4.r
                            + tw4r * t5.r
                            + tw5r * t6.r;
                        ca.i = t1.i
                            + tw1r * t2.i
                            + tw2r * t3.i
                            + tw3r * t4.i
                            + tw4r * t5.i
                            + tw5r * t6.i;
                        cb.i =
                            tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r;
                        cb.r = -(tw1i * t11.i
                            + tw2i * t10.i
                            + tw3i * t9.i
                            + tw4i * t8.i
                            + tw5i * t7.i);
                        {
                            da.r = ca.r + cb.r;
                            da.i = ca.i + cb.i;
                            db.r = ca.r - cb.r;
                            db.i = ca.i - cb.i;
                        }
                    }
                    {
                        ch[i + ido * (k + l1)].r = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.r
                            - (sign as f64) * wa[i - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1)].i = wa[i - 1 + (1 - 1) * (ido - 1)].r * da.i
                            + (sign as f64) * wa[i - 1 + (1 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 10)].r = wa[i - 1 + (10 - 1) * (ido - 1)].r * db.r
                            - (sign as f64) * wa[i - 1 + (10 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 10)].i = wa[i - 1 + (10 - 1) * (ido - 1)].r * db.i
                            + (sign as f64) * wa[i - 1 + (10 - 1) * (ido - 1)].i * db.r;
                    }
                }
                {
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                        let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                        ca.r = t1.r
                            + tw2r * t2.r
                            + tw4r * t3.r
                            + tw5r * t4.r
                            + tw3r * t5.r
                            + tw1r * t6.r;
                        ca.i = t1.i
                            + tw2r * t2.i
                            + tw4r * t3.i
                            + tw5r * t4.i
                            + tw3r * t5.i
                            + tw1r * t6.i;
                        cb.i =
                            tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r;
                        cb.r = -(tw2i * t11.i + tw4i * t10.i
                            - tw5i * t9.i
                            - tw3i * t8.i
                            - tw1i * t7.i);
                        {
                            da.r = ca.r + cb.r;
                            da.i = ca.i + cb.i;
                            db.r = ca.r - cb.r;
                            db.i = ca.i - cb.i;
                        }
                    }
                    {
                        ch[i + ido * (k + l1 * 2)].r = wa[i - 1 + (ido - 1)].r * da.r
                            - (sign as f64) * wa[i - 1 + (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1 * 2)].i = wa[i - 1 + (ido - 1)].r * da.i
                            + (sign as f64) * wa[i - 1 + (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 9)].r = wa[i - 1 + (9 - 1) * (ido - 1)].r * db.r
                            - (sign as f64) * wa[i - 1 + (9 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 9)].i = wa[i - 1 + (9 - 1) * (ido - 1)].r * db.i
                            + (sign as f64) * wa[i - 1 + (9 - 1) * (ido - 1)].i * db.r;
                    }
                }
                {
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                        let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                        ca.r = t1.r
                            + tw3r * t2.r
                            + tw5r * t3.r
                            + tw2r * t4.r
                            + tw1r * t5.r
                            + tw4r * t6.r;
                        ca.i = t1.i
                            + tw3r * t2.i
                            + tw5r * t3.i
                            + tw2r * t4.i
                            + tw1r * t5.i
                            + tw4r * t6.i;
                        cb.i =
                            tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r;
                        cb.r = -(tw3i * t11.i - tw5i * t10.i - tw2i * t9.i
                            + tw1i * t8.i
                            + tw4i * t7.i);
                        {
                            da.r = ca.r + cb.r;
                            da.i = ca.i + cb.i;
                            db.r = ca.r - cb.r;
                            db.i = ca.i - cb.i;
                        }
                    }
                    {
                        ch[i + ido * (k + l1 * 3)].r = wa[i - 1 + (3 - 1) * (ido - 1)].r * da.r
                            - (sign as f64) * wa[i - 1 + (3 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1 * 3)].i = wa[i - 1 + (3 - 1) * (ido - 1)].r * da.i
                            + (sign as f64) * wa[i - 1 + (3 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 8)].r = wa[i - 1 + (8 - 1) * (ido - 1)].r * db.r
                            - (sign as f64) * wa[i - 1 + (8 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 8)].i = wa[i - 1 + (8 - 1) * (ido - 1)].r * db.i
                            + (sign as f64) * wa[i - 1 + (8 - 1) * (ido - 1)].i * db.r;
                    }
                }
                {
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                        let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                        ca.r = t1.r
                            + tw4r * t2.r
                            + tw3r * t3.r
                            + tw1r * t4.r
                            + tw5r * t5.r
                            + tw2r * t6.r;
                        ca.i = t1.i
                            + tw4r * t2.i
                            + tw3r * t3.i
                            + tw1r * t4.i
                            + tw5r * t5.i
                            + tw2r * t6.i;
                        cb.i =
                            tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r;
                        cb.r = -(tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i
                            - tw2i * t7.i);
                        {
                            da.r = ca.r + cb.r;
                            da.i = ca.i + cb.i;
                            db.r = ca.r - cb.r;
                            db.i = ca.i - cb.i;
                        }
                    }
                    {
                        ch[i + ido * (k + l1 * 4)].r = wa[i - 1 + (4 - 1) * (ido - 1)].r * da.r
                            - (sign as f64) * wa[i - 1 + (4 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1 * 4)].i = wa[i - 1 + (4 - 1) * (ido - 1)].r * da.i
                            + (sign as f64) * wa[i - 1 + (4 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 7)].r = wa[i - 1 + (7 - 1) * (ido - 1)].r * db.r
                            - (sign as f64) * wa[i - 1 + (7 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 7)].i = wa[i - 1 + (7 - 1) * (ido - 1)].r * db.i
                            + (sign as f64) * wa[i - 1 + (7 - 1) * (ido - 1)].i * db.r;
                    }
                }
                {
                    let mut da: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut db: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        let mut ca: cmplx = cmplx { r: 0.0, i: 0.0 };
                        let mut cb: cmplx = cmplx { r: 0.0, i: 0.0 };
                        ca.r = t1.r
                            + tw5r * t2.r
                            + tw1r * t3.r
                            + tw4r * t4.r
                            + tw2r * t5.r
                            + tw3r * t6.r;
                        ca.i = t1.i
                            + tw5r * t2.i
                            + tw1r * t3.i
                            + tw4r * t4.i
                            + tw2r * t5.i
                            + tw3r * t6.i;
                        cb.i =
                            tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r;
                        cb.r = -(tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i
                            + tw3i * t7.i);
                        {
                            da.r = ca.r + cb.r;
                            da.i = ca.i + cb.i;
                            db.r = ca.r - cb.r;
                            db.i = ca.i - cb.i;
                        }
                    }
                    {
                        ch[i + ido * (k + l1 * 5)].r = wa[i - 1 + (5 - 1) * (ido - 1)].r * da.r
                            - (sign as f64) * wa[i - 1 + (5 - 1) * (ido - 1)].i * da.i;
                        ch[i + ido * (k + l1 * 5)].i = wa[i - 1 + (5 - 1) * (ido - 1)].r * da.i
                            + (sign as f64) * wa[i - 1 + (5 - 1) * (ido - 1)].i * da.r;
                    }
                    {
                        ch[i + ido * (k + l1 * 6)].r = wa[i - 1 + (6 - 1) * (ido - 1)].r * db.r
                            - (sign as f64) * wa[i - 1 + (6 - 1) * (ido - 1)].i * db.i;
                        ch[i + ido * (k + l1 * 6)].i = wa[i - 1 + (6 - 1) * (ido - 1)].r * db.i
                            + (sign as f64) * wa[i - 1 + (6 - 1) * (ido - 1)].i * db.r;
                    }
                }
            }
        }
    }
}

fn passg(
    ido: usize,
    ip: usize,
    l1: usize,
    cc: &mut [cmplx],
    ch: &mut [cmplx],
    wa: &[cmplx],
    csarr: &[cmplx],
    sign: i64,
) -> i32 {
    let cdim: usize = ip;
    let ipph: usize = (ip + 1) / 2;
    let idl1: usize = ido * l1;

    let mut wal: Vec<cmplx> = Vec::new();
    wal.push(cmplx { r: 1.0, i: 0.0 });
    for i in 1..ip {
        wal.push(cmplx {
            r: csarr[i].r,
            i: (sign as f64) * csarr[i].i,
        });
    }

    for k in 0..l1 {
        for i in 0..ido {
            ch[i + ido * k] = cc[i + ido * cdim * k];
        }
    }

    let mut j = 1;
    let mut jc = ip - 1;
    while j < ipph {
        for k in 0..l1 {
            for i in 0..ido {
                ch[i + ido * (k + l1 * j)].r =
                    cc[i + ido * (j + cdim * k)].r + cc[i + ido * (jc + cdim * k)].r;
                ch[i + ido * (k + l1 * j)].i =
                    cc[i + ido * (j + cdim * k)].i + cc[i + ido * (jc + cdim * k)].i;
                ch[i + ido * (k + l1 * jc)].r =
                    cc[i + ido * (j + cdim * k)].r - cc[i + ido * (jc + cdim * k)].r;
                ch[i + ido * (k + l1 * jc)].i =
                    cc[i + ido * (j + cdim * k)].i - cc[i + ido * (jc + cdim * k)].i;
            }
        }
        j += 1;
        jc -= 1;
    }

    for k in 0..l1 {
        for i in 0..ido {
            let mut tmp = ch[i + ido * k];
            for j in 1..ipph {
                tmp.r += ch[i + ido * (k + l1 * j)].r;
                tmp.i += ch[i + ido * (k + l1 * j)].i;
            }
            cc[i + ido * k] = tmp;
        }
    }

    let mut l = 1;
    let mut lc = ip - 1;
    while l < ipph {
        for ik in 0..idl1 {
            cc[ik + idl1 * l].r =
                ch[ik].r + wal[l].r * ch[ik + idl1].r + wal[2 * l].r * ch[ik + idl1 * 2].r;
            cc[ik + idl1 * l].i =
                ch[ik].i + wal[l].r * ch[ik + idl1].i + wal[2 * l].r * ch[ik + idl1 * 2].i;
            cc[ik + idl1 * lc].r =
                -wal[l].i * ch[ik + idl1 * (ip - 1)].i - wal[2 * l].i * ch[ik + idl1 * (ip - 2)].i;
            cc[ik + idl1 * lc].i =
                wal[l].i * ch[ik + idl1 * (ip - 1)].r + wal[2 * l].i * ch[ik + idl1 * (ip - 2)].r;
        }

        let mut iwal: usize = 2 * l;
        let mut j: usize = 3;
        let mut jc: usize = ip - 3;
        while j < ipph - 1 {
            iwal += l;
            if iwal > ip {
                iwal -= ip;
            }
            let xwal = wal[iwal];
            iwal += l;
            if iwal > ip {
                iwal -= ip;
            }
            let xwal2 = wal[iwal];
            for ik in 0..idl1 {
                cc[ik + idl1 * l].r +=
                    ch[ik + idl1 * j].r * xwal.r + ch[ik + idl1 * (j + 1)].r * xwal2.r;
                cc[ik + idl1 * l].i +=
                    ch[ik + idl1 * j].i * xwal.r + ch[ik + idl1 * (j + 1)].i * xwal2.r;
                cc[ik + idl1 * lc].r -=
                    ch[ik + idl1 * jc].i * xwal.i + ch[ik + idl1 * (jc - 1)].i * xwal2.i;
                cc[ik + idl1 * lc].i +=
                    ch[ik + idl1 * jc].r * xwal.i + ch[ik + idl1 * (jc - 1)].r * xwal2.i;
            }
            j += 2;
            jc -= 2;
        }
        while j < ipph {
            iwal += l;
            if iwal > ip {
                iwal -= ip;
            }
            let xwal = wal[iwal];
            for ik in 0..idl1 {
                cc[ik + idl1 * l].r += ch[ik + idl1 * j].r * xwal.r;
                cc[ik + idl1 * l].i += ch[ik + idl1 * j].i * xwal.r;
                cc[ik + idl1 * lc].r -= ch[ik + idl1 * jc].i * xwal.i;
                cc[ik + idl1 * lc].i += ch[ik + idl1 * jc].r * xwal.i;
            }
            j += 1;
            jc -= 1;
        }
        l += 1;
        lc -= 1;
    }

    if ido == 1 {
        let mut j = 1;
        let mut jc = ip - 1;
        while j < ipph {
            for ik in 0..idl1 {
                let t1 = cc[ik + idl1 * j];
                let t2 = cc[ik + idl1 * jc];
                {
                    cc[ik + idl1 * j].r = t1.r + t2.r;
                    cc[ik + idl1 * j].i = t1.i + t2.i;
                    cc[ik + idl1 * jc].r = t1.r - t2.r;
                    cc[ik + idl1 * jc].i = t1.i - t2.i;
                }
            }
            j += 1;
            jc -= 1;
        }
    } else {
        let mut j = 1;
        let mut jc = ip - 1;
        while j < ipph {
            for k in 0..l1 {
                let t1 = cc[ido * (k + l1 * j)];
                let t2 = cc[ido * (k + l1 * jc)];
                {
                    cc[ido * (k + l1 * j)].r = t1.r + t2.r;
                    cc[ido * (k + l1 * j)].i = t1.i + t2.i;
                    cc[ido * (k + l1 * jc)].r = t1.r - t2.r;
                    cc[ido * (k + l1 * jc)].i = t1.i - t2.i;
                }
                for i in 1..ido {
                    let mut x1: cmplx = cmplx { r: 0.0, i: 0.0 };
                    let mut x2: cmplx = cmplx { r: 0.0, i: 0.0 };
                    {
                        x1.r = cc[i + ido * (k + l1 * j)].r + cc[i + ido * (k + l1 * jc)].r;
                        x1.i = cc[i + ido * (k + l1 * j)].i + cc[i + ido * (k + l1 * jc)].i;
                        x2.r = cc[i + ido * (k + l1 * j)].r - cc[i + ido * (k + l1 * jc)].r;
                        x2.i = cc[i + ido * (k + l1 * j)].i - cc[i + ido * (k + l1 * jc)].i;
                    }
                    let mut idij = (j - 1) * (ido - 1) + i - 1;
                    {
                        cc[i + ido * (k + l1 * j)].r =
                            wa[idij].r * x1.r - (sign as f64) * wa[idij].i * x1.i;
                        cc[i + ido * (k + l1 * j)].i =
                            wa[idij].r * x1.i + (sign as f64) * wa[idij].i * x1.r;
                    }
                    idij = (jc - 1) * (ido - 1) + i - 1;
                    {
                        cc[i + ido * (k + l1 * jc)].r =
                            wa[idij].r * x2.r - (sign as f64) * wa[idij].i * x2.i;
                        cc[i + ido * (k + l1 * jc)].i =
                            wa[idij].r * x2.i + (sign as f64) * wa[idij].i * x2.r;
                    }
                }
            }
            j += 1;
            jc -= 1;
        }
    }
    return 0;
}

#[inline(never)]
fn pass_all(plan: &mut cfftp_plan_i, c: &mut [cmplx], fct: f64, sign: i32) -> i32 {
    let len = plan.length;
    if len == 1 {
        return 0;
    }
    let mut l1 = 1;
    let nf = plan.nfct;
    let mut ch: Vec<cmplx> = vec![cmplx { r: 0.0, i: 0.0 }; len];

    let c_ptr = c.as_ptr();
    let mut p1: &mut [cmplx] = &mut c[..len];
    let mut p2: &mut [cmplx] = ch.as_mut_slice();
    let mut k1 = 0;
    while k1 < nf {
        let ip: usize = plan.fct[k1].fct;
        let l2: usize = ip * l1;
        let ido: usize = len / l2;
        let wa = &plan.mem[plan.fct[k1].tw_index..];
        if ip == 4 {
            if sign > 0 {
                pass4b(ido, l1, p1, p2, wa);
            } else {
                pass4f(ido, l1, p1, p2, wa);
            };
        } else if ip == 2 {
            if sign > 0 {
                pass2b(ido, l1, p1, p2, wa);
            } else {
                pass2f(ido, l1, p1, p2, wa);
            };
        } else if ip == 3 {
            if sign > 0 {
                pass3b(ido, l1, p1, p2, wa);
            } else {
                pass3f(ido, l1, p1, p2, wa);
            };
        } else if ip == 5 {
            if sign > 0 {
                pass5b(ido, l1, p1, p2, wa);
            } else {
                pass5f(ido, l1, p1, p2, wa);
            };
        } else if ip == 7 {
            pass7(ido, l1, p1, p2, wa, sign as i64);
        } else if ip == 11 {
            pass11(ido, l1, p1, p2, wa, sign as i64);
        } else {
            let csarr = &plan.mem[plan.fct[k1].tws_index..];
            let res = passg(ido, ip, l1, p1, p2, wa, csarr, sign as i64);
            if res != 0 {
                return -1;
            }
            std::mem::swap(&mut p1, &mut p2);
        }
        std::mem::swap(&mut p1, &mut p2);
        l1 = l2;
        k1 += 1;
    }
    if p1.as_ptr() != c_ptr {
        if fct != 1.0f64 {
            for i in 0..len {
                p2[i].r = p1[i].r * fct;
                p2[i].i = p1[i].i * fct;
            }
        } else {
            p2.copy_from_slice(p1);
        }
    } else if fct != 1.0f64 {
        for i in 0..len {
            p1[i].r *= fct;
            p1[i].i *= fct;
        }
    }
    return 0;
}
#[inline(never)]
fn cfftp_forward(plan: &mut cfftp_plan_i, c: &mut [f64], fct: f64) -> i32 {
    let c = unsafe { std::mem::transmute::<&mut [f64], &mut [cmplx]>(c) };
    return pass_all(plan, c, fct, -1);
}
#[inline(never)]
fn cfftp_backward(plan: &mut cfftp_plan_i, c: &mut [f64], fct: f64) -> i32 {
    let c = unsafe { std::mem::transmute::<&mut [f64], &mut [cmplx]>(c) };
    return pass_all(plan, c, fct, 1);
}
#[inline(never)]
fn cfftp_factorize(plan: &mut cfftp_plan_i) -> i32 {
    let mut length = plan.length;
    let mut nfct: usize = 0;
    while (length % 4) == 0 {
        if nfct >= NFCT {
            return -1;
        }
        let fresh1 = nfct;
        nfct += 1;
        plan.fct[fresh1].fct = 4;
        length >>= 2;
    }
    if (length % 2) == 0 {
        length >>= 1;
        if nfct >= NFCT {
            return -1;
        }
        let fresh2 = nfct;
        nfct += 1;

        plan.fct[fresh2].fct = 2;
        let tmp = plan.fct[0].fct;
        plan.fct[0].fct = plan.fct[nfct - 1].fct;
        plan.fct[nfct - 1].fct = tmp;
    }
    let mut maxl: usize = ((length as f64).sqrt() as usize) + 1;
    let mut divisor: usize = 3;
    while length > 1 && divisor < maxl {
        if (length % divisor) == 0 {
            while (length % divisor) == 0 {
                if nfct >= NFCT {
                    return -1;
                }
                let fresh3 = nfct;
                nfct += 1;
                plan.fct[fresh3].fct = divisor;
                length /= divisor;
            }
            maxl = ((length as f64).sqrt() as usize) + 1;
        }
        divisor += 2;
    }
    if length > 1 {
        let fresh4 = nfct;
        nfct += 1;
        plan.fct[fresh4].fct = length;
    }
    plan.nfct = nfct;
    return 0;
}

#[inline(never)]
fn cfftp_twsize(plan: &mut cfftp_plan_i) -> usize {
    let mut twsize: usize = 0;
    let mut l1: usize = 1;
    let mut k: usize = 0;
    let nfct = plan.nfct;
    while k < nfct {
        let ip: usize = plan.fct[k].fct;
        let ido: usize = plan.length / (l1 * ip);
        twsize += (ip - 1) * (ido - 1);
        if ip > 11 {
            twsize += ip;
        }
        l1 *= ip;
        k += 1;
    }
    return twsize;
}

#[inline(never)]
fn cfftp_comp_twiddle(plan: &mut cfftp_plan_i) -> i32 {
    let length: usize = plan.length;
    let mut twid = vec![0.0f64; 2 * length];
    sincos_2pibyn(length, twid.as_mut_slice());
    let mut l1: usize = 1;
    let mut memofs: usize = 0;
    let mut k: usize = 0;
    let nfct: usize = plan.nfct;
    while k < nfct {
        let ip: usize = plan.fct[k].fct;
        let ido: usize = length / (l1 * ip);
        plan.fct[k].tw_index = memofs;
        memofs += (ip - 1) * (ido - 1);

        for j in 1..ip {
            let mut i: usize = 1;
            let tw = &mut plan.mem[plan.fct[k].tw_index..];
            while i < ido {
                tw[(j - 1) * (ido - 1) + i - 1].r = twid[2 * j * l1 * i];
                tw[(j - 1) * (ido - 1) + i - 1].i = twid[2 * j * l1 * i + 1];
                i += 1;
            }
        }
        if ip > 11 {
            plan.fct[k].tws_index = memofs;
            let tws = &mut plan.mem[plan.fct[k].tws_index..];
            for j in 0..ip {
                tws[j].r = twid[2 * j * l1 * ido];
                tws[j].i = twid[2 * j * l1 * ido + 1];
            }
            memofs += ip;
        }
        l1 *= ip;
        k += 1;
    }
    return 0;
}

fn make_cfftp_plan(length: usize) -> cfftp_plan {
    if length == 0 {
        return null_mut();
    }
    let tmp_fct = [cfftp_fctdata {
        fct: 0,
        tw_index: 0,
        tws_index: 0,
    }; NFCT];
    let mut tmp_cfftp_plan_i = cfftp_plan_i {
        length: 0,
        nfct: 0,
        mem: Vec::new(),
        fct: tmp_fct,
    };

    tmp_cfftp_plan_i.length = length;
    tmp_cfftp_plan_i.nfct = 0;
    let mut i: usize = 0;
    while i < NFCT {
        tmp_cfftp_plan_i.fct[i] = cfftp_fctdata {
            fct: 0,
            tw_index: 0,
            tws_index: 0,
        };
        i += 1;
    }

    if length == 1 {
        let plan: cfftp_plan = Box::into_raw(Box::new(tmp_cfftp_plan_i));
        return plan;
    }

    if cfftp_factorize(&mut tmp_cfftp_plan_i) != 0 {
        return null_mut();
    }
    let tws = cfftp_twsize(&mut tmp_cfftp_plan_i);
    tmp_cfftp_plan_i.mem = vec![cmplx { r: 0.0, i: 0.0 }; tws];

    if cfftp_comp_twiddle(&mut tmp_cfftp_plan_i) != 0 {
        return null_mut();
    }

    let plan: cfftp_plan = Box::into_raw(Box::new(tmp_cfftp_plan_i));
    return plan;
}

fn destroy_cfftp_plan(plan: cfftp_plan) {
    assert!(!plan.is_null(), "null");
    unsafe {
        let _ = Box::from_raw(plan);
    }
}

fn radf2(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 2;
    let mut k: usize = 0;
    while k < l1 {
        ch[ido * cdim * k] = cc[ido * k] + cc[ido * (k + l1)];
        ch[(ido - 1) + ido * (1 + cdim * k)] = cc[ido * k] - cc[ido * (k + l1)];
        k += 1;
    }
    if (ido & 1) == 0 {
        let mut k: usize = 0;
        while k < l1 {
            ch[ido * (1 + cdim * k)] = -cc[(ido - 1) + ido * (k + l1)];
            ch[(ido - 1) + ido * cdim * k] = cc[(ido - 1) + ido * k];
            k += 1;
        }
    }
    if ido <= 2 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic: usize = ido - i;
            let tr2: f64;
            let ti2: f64;
            tr2 = wa[i - 2] * cc[(i - 1) + ido * (k + l1)] + wa[i - 1] * cc[i + ido * (k + l1)];
            ti2 = wa[i - 2] * cc[i + ido * (k + l1)] - wa[i - 1] * cc[(i - 1) + ido * (k + l1)];
            ch[(i - 1) + ido * cdim * k] = cc[(i - 1) + ido * k] + tr2;
            ch[(ic - 1) + ido * (1 + cdim * k)] = cc[(i - 1) + ido * k] - tr2;
            ch[i + ido * cdim * k] = ti2 + cc[i + ido * k];
            ch[ic + ido * (1 + cdim * k)] = ti2 - cc[i + ido * k];
            i += 2;
        }
        k += 1;
    }
}

fn radf3(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 3;
    let taur: f64 = -0.5;
    let taui: f64 = 0.866_025_403_784_438_6;

    let mut k: usize = 0;
    while k < l1 {
        let cr2: f64 = cc[ido * (k + l1)] + cc[ido * (k + l1 * 2)];
        ch[ido * cdim * k] = cc[ido * k] + cr2;
        ch[ido * (2 + cdim * k)] = taui * (cc[ido * (k + l1 * 2)] - cc[ido * (k + l1)]);
        ch[(ido - 1) + ido * (1 + cdim * k)] = cc[ido * k] + taur * cr2;
        k += 1;
    }
    if ido == 1 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic: usize = ido - i;
            let di2: f64;
            let di3: f64;
            let dr2: f64;
            let dr3: f64;
            {
                dr2 = wa[i - 2] * cc[(i - 1) + ido * (k + l1)] + wa[i - 1] * cc[i + ido * (k + l1)];
                di2 = wa[i - 2] * cc[i + ido * (k + l1)] - wa[i - 1] * cc[(i - 1) + ido * (k + l1)];
            }
            {
                dr3 = wa[(i - 2) + (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 2)]
                    + wa[(i - 1) + (ido - 1)] * cc[i + ido * (k + l1 * 2)];
                di3 = wa[(i - 2) + (ido - 1)] * cc[i + ido * (k + l1 * 2)]
                    - wa[(i - 1) + (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 2)];
            }
            let cr2: f64 = dr2 + dr3;
            let ci2: f64 = di2 + di3;
            ch[(i - 1) + ido * cdim * k] = cc[(i - 1) + ido * k] + cr2;
            ch[i + ido * cdim * k] = cc[i + ido * k] + ci2;
            let tr2: f64 = cc[(i - 1) + ido * k] + taur * cr2;
            let ti2: f64 = cc[i + ido * k] + taur * ci2;
            let tr3: f64 = taui * (di2 - di3);
            let ti3: f64 = taui * (dr3 - dr2);
            {
                ch[(i - 1) + ido * (2 + cdim * k)] = tr2 + tr3;
                ch[(ic - 1) + ido * (1 + cdim * k)] = tr2 - tr3;
            }
            {
                ch[i + ido * (2 + cdim * k)] = ti3 + ti2;
                ch[ic + ido * (1 + cdim * k)] = ti3 - ti2;
            }
            i += 2;
        }
        k += 1;
    }
}

fn radf4(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 4;
    let hsqt2: f64 = 0.707_106_781_186_547_6;
    let mut k: usize = 0;
    while k < l1 {
        let tr1: f64;
        let tr2: f64;
        {
            tr1 = cc[ido * (k + l1 * 3)] + cc[ido * (k + l1)];
            ch[ido * (2 + cdim * k)] = cc[ido * (k + l1 * 3)] - cc[ido * (k + l1)];
        }
        {
            tr2 = cc[ido * k] + cc[ido * (k + l1 * 2)];
            ch[(ido - 1) + ido * (1 + cdim * k)] = cc[ido * k] - cc[ido * (k + l1 * 2)];
        }
        {
            ch[ido * cdim * k] = tr2 + tr1;
            ch[(ido - 1) + ido * (3 + cdim * k)] = tr2 - tr1;
        }
        k += 1;
    }
    if (ido & 1) == 0 {
        let mut k: usize = 0;
        while k < l1 {
            let ti1: f64 =
                -hsqt2 * (cc[(ido - 1) + ido * (k + l1)] + cc[(ido - 1) + ido * (k + l1 * 3)]);
            let tr1: f64 =
                hsqt2 * (cc[(ido - 1) + ido * (k + l1)] - cc[(ido - 1) + ido * (k + l1 * 3)]);
            {
                ch[(ido - 1) + ido * cdim * k] = cc[(ido - 1) + ido * k] + tr1;
                ch[(ido - 1) + ido * (2 + cdim * k)] = cc[(ido - 1) + ido * k] - tr1;
            }
            {
                ch[ido * (3 + cdim * k)] = ti1 + cc[(ido - 1) + ido * (k + l1 * 2)];
                ch[ido * (1 + cdim * k)] = ti1 - cc[(ido - 1) + ido * (k + l1 * 2)];
            }
            k += 1;
        }
    }
    if ido > 2 {
        let mut k: usize = 0;
        while k < l1 {
            let mut i: usize = 2;
            while i < ido {
                let ic: usize = ido - i;
                let cr2: f64 =
                    wa[i - 2] * cc[(i - 1) + ido * (k + l1)] + wa[i - 1] * cc[(i) + ido * (k + l1)];
                let ci2: f64 =
                    wa[i - 2] * cc[(i) + ido * (k + l1)] - wa[i - 1] * cc[(i - 1) + ido * (k + l1)];
                let cr3: f64 = wa[(i - 2) + (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 2)]
                    + wa[(i - 1) + (ido - 1)] * cc[(i) + ido * (k + l1 * 2)];
                let ci3: f64 = wa[(i - 2) + (ido - 1)] * cc[(i) + ido * (k + l1 * 2)]
                    - wa[(i - 1) + (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 2)];
                let cr4: f64 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 3)]
                    + wa[(i - 1) + (2) * (ido - 1)] * cc[(i) + ido * (k + l1 * 3)];
                let ci4: f64 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i) + ido * (k + l1 * 3)]
                    - wa[(i - 1) + (2) * (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 3)];
                let tr1: f64 = cr4 + cr2;
                let tr4: f64 = cr4 - cr2;
                let ti1: f64 = ci2 + ci4;
                let ti4: f64 = ci2 - ci4;
                let tr2: f64 = cc[(i - 1) + ido * k] + cr3;
                let tr3: f64 = cc[(i - 1) + ido * k] - cr3;
                let ti2: f64 = cc[i + ido * k] + ci3;
                let ti3: f64 = cc[i + ido * k] - ci3;
                ch[(i - 1) + ido * cdim * k] = tr2 + tr1;
                ch[(ic - 1) + ido * (3 + cdim * k)] = tr2 - tr1;
                ch[i + ido * cdim * k] = ti1 + ti2;
                ch[ic + ido * (3 + cdim * k)] = ti1 - ti2;
                ch[(i - 1) + ido * (2 + cdim * k)] = tr3 + ti4;
                ch[(ic - 1) + ido * (1 + cdim * k)] = tr3 - ti4;
                ch[i + ido * (2 + cdim * k)] = tr4 + ti3;
                ch[ic + ido * (1 + cdim * k)] = tr4 - ti3;
                i += 2;
            }
            k += 1;
        }
    }
}

fn radf5(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 5;
    let tr11: f64 = 0.309_016_994_374_947_45;
    let ti11: f64 = 0.951_056_516_295_153_5;
    let tr12: f64 = -0.809_016_994_374_947_5;
    let ti12: f64 = 0.587_785_252_292_473_1;

    let mut k: usize = 0;
    while k < l1 {
        let cr2: f64 = cc[ido * (k + l1 * 4)] + cc[ido * (k + l1)];
        let ci5: f64 = cc[ido * (k + l1 * 4)] - cc[ido * (k + l1)];
        let cr3: f64 = cc[ido * (k + l1 * 3)] + cc[ido * (k + l1 * 2)];
        let ci4: f64 = cc[ido * (k + l1 * 3)] - cc[ido * (k + l1 * 2)];
        ch[ido * cdim * k] = cc[ido * k] + cr2 + cr3;
        ch[(ido - 1) + ido * (1 + cdim * k)] = cc[ido * k] + tr11 * cr2 + tr12 * cr3;
        ch[ido * (2 + cdim * k)] = ti11 * ci5 + ti12 * ci4;
        ch[(ido - 1) + ido * (3 + cdim * k)] = cc[ido * k] + tr12 * cr2 + tr11 * cr3;
        ch[ido * (4 + cdim * k)] = ti12 * ci5 - ti11 * ci4;
        k += 1;
    }
    if ido == 1 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic: usize = ido - i;

            let dr2: f64 =
                wa[i - 2] * cc[(i - 1) + ido * (k + l1)] + wa[i - 1] * cc[i + ido * (k + l1)];
            let di2: f64 =
                wa[i - 2] * cc[i + ido * (k + l1)] - wa[i - 1] * cc[(i - 1) + ido * (k + l1)];

            let dr3: f64 = wa[(i - 2) + (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 2)]
                + wa[(i - 1) + (ido - 1)] * cc[i + ido * (k + l1 * 2)];
            let di3: f64 = wa[(i - 2) + (ido - 1)] * cc[i + ido * (k + l1 * 2)]
                - wa[(i - 1) + (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 2)];

            let dr4: f64 = wa[(i - 2) + 2 * (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 3)]
                + wa[(i - 1) + (2) * (ido - 1)] * cc[i + ido * (k + l1 * 3)];
            let di4: f64 = wa[(i - 2) + 2 * (ido - 1)] * cc[i + ido * (k + l1 * 3)]
                - wa[(i - 1) + 2 * (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 3)];

            let dr5: f64 = wa[(i - 2) + 3 * (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 4)]
                + wa[(i - 1) + 3 * (ido - 1)] * cc[i + ido * (k + l1 * 4)];
            let di5: f64 = wa[(i - 2) + 3 * (ido - 1)] * cc[(i) + ido * (k + l1 * 4)]
                - wa[(i - 1) + 3 * (ido - 1)] * cc[(i - 1) + ido * (k + l1 * 4)];

            let cr2: f64 = dr5 + dr2;
            let ci5: f64 = dr5 - dr2;

            let ci2: f64 = di2 + di5;
            let cr5: f64 = di2 - di5;

            let cr3: f64 = dr4 + dr3;
            let ci4: f64 = dr4 - dr3;

            let ci3: f64 = di3 + di4;
            let cr4: f64 = di3 - di4;

            ch[(i - 1) + ido * cdim * k] = cc[(i - 1) + ido * k] + cr2 + cr3;
            ch[i + ido * cdim * k] = cc[i + ido * k] + ci2 + ci3;
            let tr2: f64 = cc[(i - 1) + ido * k] + tr11 * cr2 + tr12 * cr3;
            let ti2: f64 = cc[i + ido * k] + tr11 * ci2 + tr12 * ci3;
            let tr3: f64 = cc[(i - 1) + ido * k] + tr12 * cr2 + tr11 * cr3;
            let ti3: f64 = cc[i + ido * k] + tr12 * ci2 + tr11 * ci3;

            let tr5: f64 = cr5 * ti11 + cr4 * ti12;
            let tr4: f64 = cr5 * ti12 - cr4 * ti11;

            let ti5: f64 = ci5 * ti11 + ci4 * ti12;
            let ti4: f64 = ci5 * ti12 - ci4 * ti11;

            ch[(i - 1) + ido * (2 + cdim * k)] = tr2 + tr5;
            ch[(ic - 1) + ido * (1 + cdim * k)] = tr2 - tr5;

            ch[i + ido * (2 + cdim * k)] = ti5 + ti2;
            ch[ic + ido * (1 + cdim * k)] = ti5 - ti2;

            ch[(i - 1) + ido * (4 + cdim * k)] = tr3 + tr4;
            ch[(ic - 1) + ido * (3 + cdim * k)] = tr3 - tr4;

            ch[i + ido * (4 + cdim * k)] = ti4 + ti3;
            ch[ic + ido * (3 + cdim * k)] = ti4 - ti3;

            i += 2;
        }
        k += 1;
    }
}

fn radfg(
    ido: usize,
    ip: usize,
    l1: usize,
    cc: &mut [f64],
    ch: &mut [f64],
    wa: &[f64],
    csarr: &[f64],
) {
    let cdim: usize = ip;
    let ipph: usize = (ip + 1) / 2;
    let idl1: usize = ido * l1;

    if ido > 1 {
        let mut j: usize = 1;
        let mut jc: usize = ip - 1;
        while j < ipph {
            let is: usize = (j - 1) * (ido - 1);
            let is2: usize = (jc - 1) * (ido - 1);
            let mut k: usize = 0;
            while k < l1 {
                let mut idij: usize = is;
                let mut idij2: usize = is2;
                let mut i: usize = 1;
                while i <= ido - 2 {
                    let t1 = cc[i + ido * (k + l1 * j)];
                    let t2 = cc[(i + 1) + ido * (k + l1 * j)];
                    let t3 = cc[i + ido * (k + l1 * jc)];
                    let t4 = cc[(i + 1) + ido * (k + l1 * jc)];
                    let x1 = wa[idij] * t1 + wa[idij + 1] * t2;
                    let x2 = wa[idij] * t2 - wa[idij + 1] * t1;
                    let x3 = wa[idij2] * t3 + wa[idij2 + 1] * t4;
                    let x4 = wa[idij2] * t4 - wa[idij2 + 1] * t3;
                    cc[i + ido * (k + l1 * j)] = x1 + x3;
                    cc[i + ido * (k + l1 * jc)] = x2 - x4;
                    cc[(i + 1) + ido * (k + l1 * j)] = x2 + x4;
                    cc[(i + 1) + ido * (k + l1 * jc)] = x3 - x1;
                    idij += 2;
                    idij2 += 2;
                    i += 2;
                }
                k += 1;
            }
            j += 1;
            jc -= 1;
        }
    }

    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let mut k: usize = 0;
        while k < l1 {
            let t1 = cc[ido * (k + l1 * j)];
            let t2 = cc[ido * (k + l1 * jc)];
            cc[ido * (k + l1 * j)] = t1 + t2;
            cc[ido * (k + l1 * jc)] = t2 - t1;
            k += 1;
        }
        j += 1;
        jc -= 1;
    }

    let mut l: usize = 1;
    let mut lc: usize = ip - 1;
    while l < ipph {
        let mut ik: usize = 0;
        while ik < idl1 {
            ch[ik + idl1 * l] =
                cc[ik] + csarr[2 * l] * cc[ik + idl1] + csarr[4 * l] * cc[ik + idl1 * 2];
            ch[ik + idl1 * lc] = csarr[2 * l + 1] * cc[ik + idl1 * (ip - 1)]
                + csarr[4 * l + 1] * cc[ik + idl1 * (ip - 2)];
            ik += 1;
        }
        let mut iang: usize = 2 * l;
        let mut j: usize = 3;
        let mut jc = ip - 3;
        while j < ipph - 3 {
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar1: f64 = csarr[2 * iang];
            let ai1: f64 = csarr[2 * iang + 1];
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar2: f64 = csarr[2 * iang];
            let ai2: f64 = csarr[2 * iang + 1];
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar3: f64 = csarr[2 * iang];
            let ai3: f64 = csarr[2 * iang + 1];
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar4: f64 = csarr[2 * iang];
            let ai4: f64 = csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                ch[ik + idl1 * l] += ar1 * cc[ik + idl1 * j]
                    + ar2 * cc[ik + idl1 * (j + 1)]
                    + ar3 * cc[ik + idl1 * (j + 2)]
                    + ar4 * cc[ik + idl1 * (j + 3)];
                ch[ik + idl1 * lc] += ai1 * cc[ik + idl1 * jc]
                    + ai2 * cc[ik + idl1 * (jc - 1)]
                    + ai3 * cc[ik + idl1 * (jc - 2)]
                    + ai4 * cc[ik + idl1 * (jc - 3)];
                ik += 1;
            }
            j += 4;
            jc -= 4;
        }
        while j < ipph - 1 {
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar1: f64 = csarr[2 * iang];
            let ai1: f64 = csarr[2 * iang + 1];
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar2: f64 = csarr[2 * iang];
            let ai2: f64 = csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                ch[ik + idl1 * l] += ar1 * cc[ik + idl1 * j] + ar2 * cc[ik + idl1 * (j + 1)];
                ch[ik + idl1 * lc] += ai1 * cc[ik + idl1 * jc] + ai2 * cc[ik + idl1 * (jc - 1)];
                ik += 1;
            }
            j += 2;
            jc -= 2;
        }
        while j < ipph {
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar: f64 = csarr[2 * iang];
            let ai: f64 = csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                ch[ik + idl1 * l] += ar * cc[ik + idl1 * j];
                ch[ik + idl1 * lc] += ai * cc[ik + idl1 * jc];
                ik += 1;
            }
            j += 1;
            jc -= 1;
        }
        l += 1;
        lc -= 1;
    }
    let mut ik: usize = 0;
    while ik < idl1 {
        ch[ik] = cc[ik];
        ik += 1;
    }

    let mut j: usize = 1;
    while j < ipph {
        let mut ik: usize = 0;
        while ik < idl1 {
            ch[ik] += cc[ik + idl1 * j];
            ik += 1;
        }
        j += 1;
    }

    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 0;
        while i < ido {
            cc[i + ido * cdim * k] = ch[i + ido * k];
            i += 1;
        }
        k += 1;
    }

    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let j2: usize = 2 * j - 1;
        let mut k: usize = 0;
        while k < l1 {
            cc[(ido - 1) + ido * (j2 + cdim * k)] = ch[ido * (k + l1 * j)];
            cc[ido * ((j2 + 1) + cdim * k)] = ch[ido * (k + l1 * jc)];
            k += 1;
        }
        j += 1;
        jc -= 1;
    }

    if ido == 1 {
        return;
    }

    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let j2: usize = 2 * j - 1;
        let mut k: usize = 0;
        while k < l1 {
            let mut i: usize = 1;
            let mut ic: usize = ido - i - 2;
            while i <= ido - 2 {
                cc[i + ido * ((j2 + 1) + cdim * k)] =
                    ch[i + ido * (k + l1 * j)] + ch[i + ido * (k + l1 * jc)];
                cc[ic + ido * (j2 + cdim * k)] =
                    ch[i + ido * (k + l1 * j)] - ch[i + ido * (k + l1 * jc)];
                cc[(i + 1) + ido * ((j2 + 1) + cdim * k)] =
                    ch[(i + 1) + ido * (k + l1 * j)] + ch[(i + 1) + ido * (k + l1 * jc)];
                cc[(ic + 1) + ido * (j2 + cdim * k)] =
                    ch[(i + 1) + ido * (k + l1 * jc)] - ch[(i + 1) + ido * (k + l1 * j)];
                i += 2;
                ic = ic.wrapping_sub(2);
            }
            k += 1;
        }
        j += 1;
        jc -= 1;
    }
}

fn radb2(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 2;

    let mut k: usize = 0;
    while k < l1 {
        ch[ido * k] = cc[ido * cdim * k] + cc[(ido - 1) + ido * (1 + cdim * k)];
        ch[ido * (k + l1)] = cc[ido * cdim * k] - cc[(ido - 1) + ido * (1 + cdim * k)];
        k += 1;
    }
    if (ido & 1) == 0 {
        let mut k: usize = 0;
        while k < l1 {
            ch[(ido - 1) + ido * k] = 2.0 * cc[(ido - 1) + ido * cdim * k];
            ch[(ido - 1) + ido * (k + l1)] = -2.0 * cc[ido * (1 + cdim * k)];
            k += 1;
        }
    }
    if ido <= 2 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic = ido - i;

            ch[(i - 1) + ido * k] =
                cc[(i - 1) + ido * cdim * k] + cc[(ic - 1) + ido * (1 + cdim * k)];
            let tr2: f64 = cc[(i - 1) + ido * cdim * k] - cc[(ic - 1) + ido * (1 + cdim * k)];

            let ti2: f64 = cc[i + ido * cdim * k] + cc[ic + ido * (1 + cdim * k)];
            ch[i + ido * k] = cc[i + ido * cdim * k] - cc[ic + ido * (1 + cdim * k)];

            ch[i + ido * (k + l1)] = wa[i - 2] * ti2 + wa[i - 1] * tr2;
            ch[(i - 1) + ido * (k + l1)] = wa[i - 2] * tr2 - wa[i - 1] * ti2;

            i += 2;
        }
        k += 1;
    }
}

fn radb3(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim = 3;
    let taur: f64 = -0.5;
    let taui: f64 = 0.866_025_403_784_438_6;

    let mut k: usize = 0;
    while k < l1 {
        let tr2: f64 = 2.0 * cc[(ido - 1) + ido * (1 + cdim * k)];
        let cr2: f64 = cc[ido * cdim * k] + taur * tr2;
        ch[ido * k] = cc[ido * cdim * k] + tr2;
        let ci3: f64 = 2.0 * taui * cc[ido * (2 + cdim * k)];

        ch[ido * (k + l1 * 2)] = cr2 + ci3;
        ch[ido * (k + l1)] = cr2 - ci3;

        k += 1;
    }
    if ido == 1 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic = ido - i;
            let tr2: f64 = cc[(i - 1) + ido * (2 + cdim * k)] + cc[(ic - 1) + ido * (1 + cdim * k)];
            let ti2: f64 = cc[i + ido * (2 + cdim * k)] - cc[ic + ido * (1 + cdim * k)];
            let cr2: f64 = cc[(i - 1) + ido * cdim * (k)] + taur * tr2;
            let ci2: f64 = cc[i + ido * cdim * k] + taur * ti2;
            ch[(i - 1) + ido * k] = cc[(i - 1) + ido * cdim * k] + tr2;
            ch[i + ido * k] = cc[i + ido * cdim * k] + ti2;
            let cr3: f64 =
                taui * (cc[(i - 1) + ido * (2 + cdim * k)] - cc[(ic - 1) + ido * (1 + cdim * k)]);
            let ci3: f64 = taui * (cc[i + ido * (2 + cdim * k)] + cc[ic + ido * (1 + cdim * k)]);

            let dr3: f64 = cr2 + ci3;
            let dr2: f64 = cr2 - ci3;

            let di2: f64 = ci2 + cr3;
            let di3: f64 = ci2 - cr3;

            ch[i + ido * (k + l1)] = wa[i - 2] * di2 + wa[i - 1] * dr2;
            ch[(i - 1) + ido * ((k) + l1)] = wa[i - 2] * dr2 - wa[i - 1] * di2;

            ch[i + ido * (k + l1 * 2)] =
                wa[(i - 2) + (ido - 1)] * di3 + wa[(i - 1) + (ido - 1)] * dr3;
            ch[(i - 1) + ido * (k + l1 * 2)] =
                wa[(i - 2) + (ido - 1)] * dr3 - wa[(i - 1) + (ido - 1)] * di3;

            i += 2;
        }
        k += 1;
    }
}

fn radb4(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 4;
    let sqrt2 = 1.414_213_562_373_095_1;

    let mut k: usize = 0;
    while k < l1 {
        let tr2: f64 = cc[ido * cdim * k] + cc[(ido - 1) + ido * (3 + cdim * k)];
        let tr1: f64 = cc[ido * cdim * k] - cc[(ido - 1) + ido * (3 + cdim * k)];

        let tr3: f64 = 2.0 * cc[(ido - 1) + ido * (1 + cdim * k)];
        let tr4: f64 = 2.0 * cc[ido * (2 + cdim * k)];
        {
            ch[ido * k] = tr2 + tr3;
            ch[ido * (k + l1 * 2)] = tr2 - tr3;
        }

        ch[ido * (k + l1 * 3)] = tr1 + tr4;
        ch[ido * (k + l1)] = tr1 - tr4;

        k += 1;
    }
    if (ido & 1) == 0 {
        let mut k: usize = 0;
        while k < l1 {
            let ti1: f64 = cc[ido * (3 + cdim * k)] + cc[ido * (1 + cdim * k)];
            let ti2: f64 = cc[ido * (3 + cdim * k)] - cc[ido * (1 + cdim * k)];

            let tr2: f64 = cc[(ido - 1) + ido * cdim * k] + cc[(ido - 1) + ido * (2 + cdim * k)];
            let tr1: f64 = cc[(ido - 1) + ido * cdim * k] - cc[(ido - 1) + ido * (2 + cdim * k)];

            ch[(ido - 1) + ido * k] = tr2 + tr2;
            ch[(ido - 1) + ido * (k + l1)] = sqrt2 * (tr1 - ti1);
            ch[(ido - 1) + ido * (k + l1 * 2)] = ti2 + ti2;
            ch[(ido - 1) + ido * (k + l1 * 3)] = -sqrt2 * (tr1 + ti1);
            k += 1;
        }
    }
    if ido <= 2 {
        return;
    };
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic: usize = ido - i;

            let tr2: f64 = cc[(i - 1) + ido * cdim * k] + cc[(ic - 1) + ido * (3 + cdim * k)];
            let tr1: f64 = cc[(i - 1) + ido * cdim * k] - cc[(ic - 1) + ido * (3 + cdim * k)];

            let ti1: f64 = cc[i + ido * cdim * k] + cc[ic + ido * (3 + cdim * k)];
            let ti2: f64 = cc[i + ido * cdim * k] - cc[ic + ido * (3 + cdim * k)];

            let tr4: f64 = cc[i + ido * (2 + cdim * k)] + cc[ic + ido * (1 + cdim * k)];
            let ti3: f64 = cc[i + ido * (2 + cdim * k)] - cc[ic + ido * (1 + cdim * k)];

            let tr3: f64 = cc[(i - 1) + ido * (2 + cdim * k)] + cc[(ic - 1) + ido * (1 + cdim * k)];
            let ti4: f64 = cc[(i - 1) + ido * (2 + cdim * k)] - cc[(ic - 1) + ido * (1 + cdim * k)];

            ch[(i - 1) + ido * k] = tr2 + tr3;
            let cr3: f64 = tr2 - tr3;

            ch[i + ido * k] = ti2 + ti3;
            let ci3: f64 = ti2 - ti3;

            let cr4: f64 = tr1 + tr4;
            let cr2: f64 = tr1 - tr4;

            let ci2: f64 = ti1 + ti4;
            let ci4: f64 = ti1 - ti4;

            ch[i + ido * (k + l1)] = wa[i - 2] * ci2 + wa[i - 1] * cr2;
            ch[(i - 1) + ido * (k + l1)] = wa[i - 2] * cr2 - wa[i - 1] * ci2;

            ch[i + ido * (k + l1 * 2)] =
                wa[(i - 2) + (ido - 1)] * ci3 + wa[(i - 1) + (ido - 1)] * cr3;
            ch[(i - 1) + ido * (k + l1 * 2)] =
                wa[(i - 2) + (ido - 1)] * cr3 - wa[(i - 1) + (ido - 1)] * ci3;

            ch[i + ido * (k + l1 * 3)] =
                wa[(i - 2) + 2 * (ido - 1)] * ci4 + wa[(i - 1) + 2 * (ido - 1)] * cr4;
            ch[(i - 1) + ido * (k + l1 * 3)] =
                wa[(i - 2) + 2 * (ido - 1)] * cr4 - wa[(i - 1) + 2 * (ido - 1)] * ci4;

            i += 2;
        }
        k += 1;
    }
}

fn radb5(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 5;
    let tr11: f64 = 0.309_016_994_374_947_45;
    let ti11: f64 = 0.951_056_516_295_153_5;
    let tr12: f64 = -0.809_016_994_374_947_5;
    let ti12: f64 = 0.587_785_252_292_473_1;

    let mut k: usize = 0;
    while k < l1 {
        let ti5 = cc[ido * (2 + cdim * k)] + cc[ido * (2 + cdim * k)];
        let ti4 = cc[ido * (4 + cdim * k)] + cc[ido * (4 + cdim * k)];
        let tr2 = cc[(ido - 1) + ido * (1 + cdim * k)] + cc[(ido - 1) + ido * (1 + cdim * k)];
        let tr3 = cc[(ido - 1) + ido * (3 + cdim * k)] + cc[(ido - 1) + ido * (3 + cdim * k)];
        ch[ido * k] = cc[ido * cdim * k] + tr2 + tr3;
        let cr2 = cc[ido * cdim * k] + tr11 * tr2 + tr12 * tr3;
        let cr3 = cc[ido * cdim * k] + tr12 * tr2 + tr11 * tr3;

        let ci5: f64 = ti5 * ti11 + ti4 * ti12;
        let ci4: f64 = ti5 * ti12 - ti4 * ti11;

        ch[ido * (k + l1 * 4)] = cr2 + ci5;
        ch[ido * (k + l1)] = cr2 - ci5;

        ch[ido * (k + l1 * 3)] = cr3 + ci4;
        ch[ido * (k + l1 * 2)] = cr3 - ci4;

        k += 1;
    }
    if ido == 1 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic: usize = ido - i;
            let tr2: f64 = cc[(i - 1) + ido * (2 + cdim * k)] + cc[(ic - 1) + ido * (1 + cdim * k)];
            let tr5: f64 = cc[(i - 1) + ido * (2 + cdim * k)] - cc[(ic - 1) + ido * (1 + cdim * k)];

            let ti5: f64 = cc[i + ido * (2 + cdim * k)] + cc[ic + ido * (1 + cdim * k)];
            let ti2: f64 = cc[i + ido * (2 + cdim * k)] - cc[ic + ido * (1 + cdim * k)];

            let tr3: f64 = cc[(i - 1) + ido * (4 + cdim * k)] + cc[(ic - 1) + ido * (3 + cdim * k)];
            let tr4: f64 = cc[(i - 1) + ido * (4 + cdim * k)] - cc[(ic - 1) + ido * (3 + cdim * k)];

            let ti4: f64 = cc[i + ido * (4 + cdim * k)] + cc[ic + ido * (3 + cdim * k)];
            let ti3: f64 = cc[i + ido * (4 + cdim * k)] - cc[ic + ido * (3 + cdim * k)];

            ch[(i - 1) + ido * k] = cc[(i - 1) + ido * cdim * k] + tr2 + tr3;
            ch[i + ido * k] = cc[i + ido * cdim * k] + ti2 + ti3;
            let cr2: f64 = cc[(i - 1) + ido * cdim * k] + tr11 * tr2 + tr12 * tr3;
            let ci2: f64 = cc[i + ido * cdim * k] + tr11 * ti2 + tr12 * ti3;
            let cr3: f64 = cc[(i - 1) + ido * cdim * k] + tr12 * tr2 + tr11 * tr3;
            let ci3: f64 = cc[i + ido * cdim * k] + tr12 * ti2 + tr11 * ti3;

            let cr5: f64 = tr5 * ti11 + tr4 * ti12;
            let cr4: f64 = tr5 * ti12 - tr4 * ti11;

            let ci5: f64 = ti5 * ti11 + ti4 * ti12;
            let ci4: f64 = ti5 * ti12 - ti4 * ti11;

            let dr4: f64 = cr3 + ci4;
            let dr3: f64 = cr3 - ci4;

            let di3: f64 = ci3 + cr4;
            let di4: f64 = ci3 - cr4;

            let dr5: f64 = cr2 + ci5;
            let dr2: f64 = cr2 - ci5;

            let di2: f64 = ci2 + cr5;
            let di5: f64 = ci2 - cr5;

            ch[i + ido * (k + l1)] = wa[i - 2] * di2 + wa[i - 1] * dr2;
            ch[(i - 1) + ido * (k + l1)] = wa[i - 2] * dr2 - wa[i - 1] * di2;

            ch[i + ido * (k + l1 * 2)] =
                wa[(i - 2) + (ido - 1)] * di3 + wa[(i - 1) + (ido - 1)] * dr3;
            ch[(i - 1) + ido * (k + l1 * 2)] =
                wa[(i - 2) + (ido - 1)] * dr3 - wa[(i - 1) + (ido - 1)] * di3;

            ch[i + ido * (k + l1 * 3)] =
                wa[(i - 2) + 2 * (ido - 1)] * di4 + wa[(i - 1) + 2 * (ido - 1)] * dr4;
            ch[(i - 1) + ido * (k + l1 * 3)] =
                wa[(i - 2) + 2 * (ido - 1)] * dr4 - wa[(i - 1) + 2 * (ido - 1)] * di4;

            ch[i + ido * (k + l1 * 4)] =
                wa[(i - 2) + 3 * (ido - 1)] * di5 + wa[(i - 1) + 3 * (ido - 1)] * dr5;
            ch[(i - 1) + ido * (k + l1 * 4)] =
                wa[(i - 2) + 3 * (ido - 1)] * dr5 - wa[(i - 1) + 3 * (ido - 1)] * di5;

            i += 2;
        }
        k += 1;
    }
}

fn radbg(
    ido: usize,
    ip: usize,
    l1: usize,
    cc: &mut [f64],
    ch: &mut [f64],
    wa: &[f64],
    csarr: &[f64],
) {
    let cdim: usize = ip;
    let ipph: usize = (ip + 1) / 2;
    let idl1: usize = ido * l1;

    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 0;
        while i < ido {
            ch[i + ido * k] = cc[i + ido * cdim * k];
            i += 1;
        }
        k += 1;
    }
    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let j2: usize = 2 * j - 1;
        let mut k: usize = 0;
        while k < l1 {
            ch[ido * (k + l1 * j)] = 2.0 * cc[(ido - 1) + ido * (j2 + cdim * k)];
            ch[ido * (k + l1 * jc)] = 2.0 * cc[ido * ((j2 + 1) + cdim * k)];
            k += 1;
        }
        j += 1;
        jc -= 1;
    }

    if ido != 1 {
        let mut j: usize = 1;
        let mut jc: usize = ip - 1;
        while j < ipph {
            let j2: usize = 2 * j - 1;
            let mut k: usize = 0;
            while k < l1 {
                let mut i: usize = 1;
                let mut ic: usize = ido - i - 2;
                while i <= ido - 2 {
                    ch[i + ido * (k + l1 * j)] =
                        cc[i + ido * ((j2 + 1) + cdim * k)] + cc[ic + ido * (j2 + cdim * k)];
                    ch[i + ido * (k + l1 * jc)] =
                        cc[i + ido * ((j2 + 1) + cdim * k)] - cc[ic + ido * (j2 + cdim * k)];
                    ch[(i + 1) + ido * (k + l1 * j)] = cc[(i + 1) + ido * ((j2 + 1) + cdim * k)]
                        - cc[(ic + 1) + ido * (j2 + cdim * k)];
                    ch[(i + 1) + ido * (k + l1 * jc)] = cc[(i + 1) + ido * ((j2 + 1) + cdim * k)]
                        + cc[(ic + 1) + ido * (j2 + cdim * k)];
                    i += 2;
                    ic = ic.wrapping_sub(2);
                }
                k += 1;
            }
            j += 1;
            jc -= 1;
        }
    }
    let mut l: usize = 1;
    let mut lc = ip - 1;
    while l < ipph {
        let mut ik: usize = 0;
        while ik < idl1 {
            cc[ik + idl1 * l] =
                ch[ik] + csarr[2 * l] * ch[ik + idl1] + csarr[4 * l] * ch[ik + idl1 * 2];
            cc[ik + idl1 * lc] = csarr[2 * l + 1] * ch[ik + idl1 * (ip - 1)]
                + csarr[4 * l + 1] * ch[ik + idl1 * (ip - 2)];
            ik += 1;
        }
        let mut iang: usize = 2 * l;
        let mut j: usize = 3;
        let mut jc = ip - 3;
        while j < ipph - 3 {
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar1 = csarr[2 * iang];
            let ai1 = csarr[2 * iang + 1];
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar2 = csarr[2 * iang];
            let ai2 = csarr[2 * iang + 1];
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar3 = csarr[2 * iang];
            let ai3 = csarr[2 * iang + 1];
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar4 = csarr[2 * iang];
            let ai4 = csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                cc[ik + idl1 * l] += ar1 * ch[ik + idl1 * j]
                    + ar2 * ch[ik + idl1 * (j + 1)]
                    + ar3 * ch[ik + idl1 * (j + 2)]
                    + ar4 * ch[ik + idl1 * (j + 3)];
                cc[ik + idl1 * lc] += ai1 * ch[ik + idl1 * jc]
                    + ai2 * ch[ik + idl1 * (jc - 1)]
                    + ai3 * ch[ik + idl1 * (jc - 2)]
                    + ai4 * ch[ik + idl1 * (jc - 3)];
                ik += 1;
            }
            j += 4;
            jc -= 4;
        }
        while j < ipph - 1 {
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar1 = csarr[2 * iang];
            let ai1 = csarr[2 * iang + 1];
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar2 = csarr[2 * iang];
            let ai2 = csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                cc[ik + idl1 * l] += ar1 * ch[ik + idl1 * j] + ar2 * ch[ik + idl1 * (j + 1)];
                cc[ik + idl1 * lc] += ai1 * ch[ik + idl1 * jc] + ai2 * ch[ik + idl1 * (jc - 1)];
                ik += 1;
            }
            j += 2;
            jc -= 2;
        }
        while j < ipph {
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let war = csarr[2 * iang];
            let wai = csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                cc[ik + idl1 * l] += war * ch[ik + idl1 * j];
                cc[ik + idl1 * lc] += wai * ch[ik + idl1 * jc];
                ik += 1;
            }
            j += 1;
            jc -= 1;
        }
        l += 1;
        lc -= 1;
    }
    let mut j: usize = 1;
    while j < ipph {
        let mut ik: usize = 0;
        while ik < idl1 {
            ch[ik] += ch[ik + idl1 * j];
            ik += 1;
        }
        j += 1;
    }
    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let mut k: usize = 0;
        while k < l1 {
            ch[ido * (k + l1 * j)] = cc[ido * (k + l1 * j)] - cc[ido * (k + l1 * jc)];
            ch[ido * (k + l1 * jc)] = cc[ido * (k + l1 * j)] + cc[ido * (k + l1 * jc)];
            k += 1;
        }
        j += 1;
        jc -= 1;
    }
    if ido == 1 {
        return;
    }

    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let mut k: usize = 0;
        while k < l1 {
            let mut i: usize = 1;
            while i <= ido - 2 {
                ch[i + ido * (k + l1 * j)] =
                    cc[i + ido * (k + l1 * j)] - cc[(i + 1) + ido * (k + l1 * jc)];
                ch[i + ido * (k + l1 * jc)] =
                    cc[i + ido * (k + l1 * j)] + cc[(i + 1) + ido * (k + l1 * jc)];
                ch[(i + 1) + ido * (k + l1 * j)] =
                    cc[(i + 1) + ido * (k + l1 * j)] + cc[i + ido * (k + l1 * jc)];
                ch[(i + 1) + ido * (k + l1 * jc)] =
                    cc[(i + 1) + ido * (k + l1 * j)] - cc[i + ido * (k + l1 * jc)];
                i += 2;
            }
            k += 1;
        }
        j += 1;
        jc -= 1;
    }

    let mut j: usize = 1;
    while j < ip {
        let is: usize = (j - 1) * (ido - 1);
        let mut k: usize = 0;
        while k < l1 {
            let mut idij: usize = is;
            let mut i: usize = 1;
            while i <= ido - 2 {
                let t1 = ch[i + ido * (k + l1 * j)];
                let t2 = ch[(i + 1) + ido * (k + l1 * j)];
                ch[(i) + ido * (k + l1 * j)] = wa[idij] * t1 - wa[idij + 1] * t2;
                ch[(i + 1) + ido * (k + l1 * j)] = wa[idij] * t2 + wa[idij + 1] * t1;
                idij += 2;
                i += 2;
            }
            k += 1;
        }
        j += 1;
    }
}

fn copy_and_norm(c: &mut [f64], p1: &mut [f64], n: usize, fct: f64) {
    if p1 != c {
        if fct != 1.0 {
            let mut i: usize = 0;
            while i < n {
                c[i] = p1[i] * fct;
                i += 1;
            }
        } else {
            c[..n].copy_from_slice(&p1[..n]);
        }
    } else if fct != 1.0 {
        let mut i: usize = 0;
        while i < n {
            c[i] *= fct;
            i += 1;
        }
    }
}

fn rfftp_forward(plan: &mut rfftp_plan_i, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = plan.length;
    if n == 1 {
        return 0;
    }
    let mut l1: usize = n;
    let nf: usize = plan.nfct;
    let mut ch = vec![0.0f64; n];
    let mut p1: &mut [f64] = &mut c[..n];
    let mut p2: &mut [f64] = ch.as_mut_slice();
    let mut k1: usize = 0;
    while k1 < nf {
        let k: usize = nf - k1 - 1;
        let ip: usize = plan.fct[k].fct;
        let ido: usize = n / l1;
        let wa = &plan.mem[plan.fct[k].tw_index..];
        l1 /= ip;
        if ip == 4 {
            radf4(ido, l1, p1, p2, wa);
        } else if ip == 2 {
            radf2(ido, l1, p1, p2, wa);
        } else if ip == 3 {
            radf3(ido, l1, p1, p2, wa);
        } else if ip == 5 {
            radf5(ido, l1, p1, p2, wa);
        } else {
            let csarr = &plan.mem[plan.fct[k].tws_index..];
            radfg(ido, ip, l1, p1, p2, wa, csarr);
            std::mem::swap(&mut p1, &mut p2);
        }
        std::mem::swap(&mut p1, &mut p2);
        k1 += 1;
    }
    copy_and_norm(p2, p1, n, fct);
    return 0;
}

fn rfftp_backward(plan: &mut rfftp_plan_i, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = plan.length;
    if n == 1 {
        return 0;
    }
    let mut l1: usize = 1;
    let nf: usize = plan.nfct;
    let mut ch = vec![0.0f64; n];
    let c_ptr = c.as_mut_ptr();
    let mut p1: &mut [f64] = &mut c[..n];
    let mut p2: &mut [f64] = ch.as_mut_slice();
    let mut k: usize = 0;
    while k < nf {
        let ip: usize = plan.fct[k].fct;
        let ido: usize = n / (ip * l1);
        let tw = plan.mem.as_slice();
        let wa = &tw[plan.fct[k].tw_index..];
        if ip == 4 {
            radb4(ido, l1, p1, p2, wa);
        } else if ip == 2 {
            radb2(ido, l1, p1, p2, wa);
        } else if ip == 3 {
            radb3(ido, l1, p1, p2, wa);
        } else if ip == 5 {
            radb5(ido, l1, p1, p2, wa);
        } else {
            let tws = plan.mem.as_slice();
            let csarr = &tws[plan.fct[k].tws_index..];
            radbg(ido, ip, l1, p1, p2, wa, csarr);
        }
        std::mem::swap(&mut p1, &mut p2);
        l1 *= ip;
        k += 1;
    }
    unsafe {
        copy_and_norm(from_raw_parts_mut(c_ptr, n), p1, n, fct);
    }
    return 0;
}

fn rfftp_factorize(plan: &mut rfftp_plan_i) -> i32 {
    let mut length = plan.length;
    let mut nfct: usize = 0;
    while (length % 4) == 0 {
        if nfct >= NFCT {
            return -1;
        }
        let fresh5 = nfct;
        nfct += 1;
        plan.fct[fresh5].fct = 4;
        length >>= 2;
    }
    if (length % 2) == 0 {
        length >>= 1;
        if nfct >= NFCT {
            return -1;
        }
        let fresh6 = nfct;
        nfct += 1;
        plan.fct[fresh6].fct = 2;
        let tmp: usize = plan.fct[0].fct;
        plan.fct[0].fct = plan.fct[nfct - 1].fct;
        plan.fct[nfct - 1].fct = tmp;
    }
    let mut maxl: usize = ((length as f64).sqrt() as usize) + 1;
    let mut divisor: usize = 3;
    while length > 1 && divisor < maxl {
        if (length % divisor) == 0 {
            while (length % divisor) == 0 {
                if nfct >= NFCT {
                    return -1;
                }
                let fresh7 = nfct;
                nfct += 1;
                plan.fct[fresh7].fct = divisor;
                length /= divisor;
            }
            maxl = ((length as f64).sqrt() as usize) + 1;
        }
        divisor += 2;
    }
    if length > 1 {
        let fresh8 = nfct;
        nfct += 1;
        plan.fct[fresh8].fct = length;
    }
    plan.nfct = nfct;
    return 0;
}

fn rfftp_twsize(plan: &mut rfftp_plan_i) -> usize {
    let mut twsize: usize = 0;
    let mut l1: usize = 1;
    for k in 0..plan.nfct {
        let ip = plan.fct[k].fct;
        let ido = plan.length / (l1 * ip);
        twsize += (ip - 1) * (ido - 1);
        if ip > 5 {
            twsize += 2 * ip;
        }
        l1 *= ip;
    }
    return twsize;
}

#[inline(never)]
fn rfftp_comp_twiddle(plan: &mut rfftp_plan_i) -> i32 {
    let length: usize = plan.length;
    let mut twid = vec![0.0f64; 2 * length];
    sincos_2pibyn_half(length, twid.as_mut_slice());
    let mut l1: usize = 1;
    let mut k: usize = 0;
    let nfct = plan.nfct;
    let mut index = 0;
    while k < nfct {
        let ip = plan.fct[k].fct;
        let ido: usize = length / (l1 * ip);
        if k < nfct - 1 {
            plan.fct[k].tw_index = index;
            let mut j: usize = 1;
            while j < ip {
                let mut i: usize = 1;
                let tw = &mut plan.mem[plan.fct[k].tw_index..];
                while i <= ((ido - 1) / 2) {
                    tw[(j - 1) * (ido - 1) + 2 * i - 2] = twid[2 * j * l1 * i];
                    tw[(j - 1) * (ido - 1) + 2 * i - 1] = twid[2 * j * l1 * i + 1];
                    i += 1;
                }
                j += 1;
            }
            index += (ip - 1) * (ido - 1);
        }
        if ip > 5 {
            plan.fct[k].tws_index = index;
            let tws = &mut plan.mem[index..];
            tws[0] = 1.0f64;
            tws[1] = 0.0f64;
            let mut i: usize = 1;
            while i <= (ip >> 1) {
                tws[2 * i] = twid[2 * i * (length / ip)];
                tws[2 * i + 1] = twid[2 * i * (length / ip) + 1];
                tws[2 * (ip - i)] = twid[2 * i * (length / ip)];
                tws[2 * (ip - i) + 1] = -twid[2 * i * (length / ip) + 1];
                i += 1;
            }
            index += 2 * ip;
        }
        l1 *= ip;
        k += 1;
    }
    return 0;
}

fn make_rfftp_plan(length: usize) -> rfftp_plan {
    if length == 0 {
        return null_mut();
    }
    let tmp_fct = [rfftp_fctdata {
        fct: 0,
        tw_index: 0,
        tws_index: 0,
    }; NFCT];
    let mut tmp_rfftp_plan_i = rfftp_plan_i {
        length,
        nfct: 0,
        mem: Vec::new(),
        fct: tmp_fct,
    };
    tmp_rfftp_plan_i.length = length;
    tmp_rfftp_plan_i.nfct = 0;

    for i in 0..NFCT {
        let init = rfftp_fctdata {
            fct: 0,
            tw_index: 0,
            tws_index: 0,
        };

        tmp_rfftp_plan_i.fct[i] = init;
    }
    if length == 1 {
        let plan: rfftp_plan = Box::into_raw(Box::new(tmp_rfftp_plan_i));
        return plan;
    }

    if rfftp_factorize(&mut tmp_rfftp_plan_i) != 0 {
        return null_mut();
    }
    let tws = rfftp_twsize(&mut tmp_rfftp_plan_i);
    tmp_rfftp_plan_i.mem = vec![0.0f64; tws];
    if rfftp_comp_twiddle(&mut tmp_rfftp_plan_i) != 0 {
        return null_mut();
    }

    let plan: rfftp_plan = Box::into_raw(Box::new(tmp_rfftp_plan_i));
    return plan;
}

fn destroy_rfftp_plan(plan: rfftp_plan) {
    unsafe {
        let _ = Box::from_raw(plan);
    }
}

fn make_fftblue_plan(length: usize) -> fftblue_plan {
    let mut tmp_fftblue_plan_i = fftblue_plan_i {
        n: 0,
        n2: 0,
        plan: null_mut(),
        mem: Vec::new(),
    };

    let n2 = good_size(length * 2 - 1);

    tmp_fftblue_plan_i.n = length;
    tmp_fftblue_plan_i.n2 = n2;

    tmp_fftblue_plan_i.mem = vec![0.0f64; 2 * length + 2 * n2];

    let mut tmp = vec![0.0f64; 4 * length];
    sincos_2pibyn(2 * length, tmp.as_mut_slice());
    let bk = tmp_fftblue_plan_i.mem.as_mut_slice();
    bk[0] = 1.0;
    bk[1] = 0.0;
    let mut coeff: usize = 0;
    let mut m: usize = 1;
    while m < length {
        coeff = coeff + m * 2 - 1;
        if coeff >= length * 2 {
            coeff -= length * 2;
        }
        bk[m * 2] = tmp[2 * coeff];
        bk[m * 2 + 1] = tmp[2 * coeff + 1];
        m += 1;
    }
    let xn2: f64 = 1.0f64 / n2 as f64;
    let bkf_index = length * 2;
    bk[bkf_index] = bk[0] * xn2;
    bk[bkf_index + 1] = bk[1] * xn2;
    m = 2;
    while m < length * 2 {
        bk[bkf_index + (2 * n2 - m)] = bk[m] * xn2;
        bk[bkf_index + m] = bk[bkf_index + (2 * n2 - m)];

        bk[bkf_index + (2 * n2 - m + 1)] = bk[m + 1] * xn2;
        bk[bkf_index + m + 1] = bk[bkf_index + (2 * n2 - m + 1)];
        m += 2;
    }

    m = length * 2;

    while m <= (2 * n2 - 2 * (length) + 1) {
        bk[bkf_index + m] = 0.0f64;
        m += 1;
    }
    tmp_fftblue_plan_i.plan = make_cfftp_plan(n2);
    if tmp_fftblue_plan_i.plan.is_null() {
        return null_mut();
    }
    let cfftp_plan_ref = unsafe { &mut *tmp_fftblue_plan_i.plan };
    if cfftp_forward(cfftp_plan_ref, &mut bk[bkf_index..], 1.0f64) != 0 {
        return null_mut();
    }

    let plan: fftblue_plan = Box::into_raw(Box::new(tmp_fftblue_plan_i));
    return plan;
}

fn destroy_fftblue_plan(plan: fftblue_plan) {
    assert!(!plan.is_null(), "null");
    unsafe {
        destroy_cfftp_plan((*plan).plan);
        let _ = Box::from_raw(plan);
    }
}

fn fftblue_fft(plan: &mut fftblue_plan_i, c: &mut [f64], isign: i32, fct: f64) -> i32 {
    let n: usize = plan.n;
    let n2: usize = plan.n2;
    let mut akf = vec![0.0f64; 2 * n2];
    let bk = plan.mem.as_slice();
    let bkf_index = n * 2;
    if isign > 0 {
        let mut m: usize = 0;
        while m < (n * 2) {
            akf[m] = c[m] * bk[m] - c[m + 1] * bk[m + 1];
            akf[m + 1] = c[m] * bk[m + 1] + c[m + 1] * bk[m];
            m += 2;
        }
    } else {
        let mut m: usize = 0;
        while m < (n * 2) {
            akf[m] = c[m] * bk[m] + c[m + 1] * bk[m + 1];
            akf[m + 1] = -c[m] * bk[m + 1] + c[m + 1] * bk[m];
            m += 2;
        }
    }
    let mut m: usize = n * 2;
    while m < (n2 * 2) {
        akf[m] = 0.0;
        m += 1;
    }
    let cfftp_plan_ref = unsafe { &mut *plan.plan };
    let res = cfftp_forward(cfftp_plan_ref, akf.as_mut_slice(), fct);
    if res != 0 {
        return -1;
    }
    if isign > 0 {
        let mut m: usize = 0;
        while m < (n2 * 2) {
            let im: f64 = -akf[m] * bk[bkf_index + m + 1] + akf[m + 1] * bk[bkf_index + m];
            akf[m] = akf[m] * bk[bkf_index + m] + akf[m + 1] * bk[bkf_index + m + 1];
            akf[m + 1] = im;
            m += 2;
        }
    } else {
        let mut m: usize = 0;
        while m < (n2 * 2) {
            let im: f64 = akf[m] * bk[bkf_index + m + 1] + akf[m + 1] * bk[bkf_index + m];
            akf[m] = akf[m] * bk[bkf_index + m] - akf[m + 1] * bk[bkf_index + m + 1];
            akf[m + 1] = im;
            m += 2;
        }
    }
    let res = cfftp_backward(cfftp_plan_ref, akf.as_mut_slice(), 1.0f64);
    if res != 0 {
        return -1;
    }
    if isign > 0 {
        let mut m: usize = 0;
        while m < (n * 2) {
            c[m] = bk[m] * akf[m] - bk[m + 1] * akf[m + 1];
            c[m + 1] = bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
            m += 2;
        }
    } else {
        let mut m: usize = 0;
        while m < (n * 2) {
            c[m] = bk[m] * akf[m] + bk[m + 1] * akf[m + 1];
            c[m + 1] = -bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
            m += 2;
        }
    }
    return 0;
}

fn cfftblue_backward(plan: &mut fftblue_plan_i, c: &mut [f64], fct: f64) -> i32 {
    return fftblue_fft(plan, c, 1, fct);
}

fn cfftblue_forward(plan: &mut fftblue_plan_i, c: &mut [f64], fct: f64) -> i32 {
    return fftblue_fft(plan, c, -1, fct);
}

fn rfftblue_backward(plan: &mut fftblue_plan_i, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = plan.n;
    let mut tmp = vec![0.0f64; 2 * n];
    tmp[0] = c[0];
    tmp[1] = 0.0f64;
    tmp[2..n + 1].copy_from_slice(&c[1..]);
    if n & 1 == 0 {
        tmp[n + 1] = 0.0f64;
    }
    let mut m: usize = 2;
    while m < n {
        tmp[2 * n - m] = tmp[m];
        tmp[2 * n - m + 1] = -tmp[m + 1];
        m += 2;
    }
    let res = fftblue_fft(plan, tmp.as_mut_slice(), 1, fct);
    if res != 0 {
        return -1;
    }

    for m in 0..n {
        c[m] = tmp[2 * m];
    }
    return 0;
}

fn rfftblue_forward(plan: &mut fftblue_plan_i, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = plan.n;
    let mut tmp = vec![0.0f64; 2 * n];
    let mut m: usize = 0;
    while m < n {
        tmp[2 * m] = c[m];
        tmp[2 * m + 1] = 0.0f64;
        m += 1;
    }
    let res = fftblue_fft(plan, tmp.as_mut_slice(), -1, fct);
    if res != 0 {
        return -1;
    }
    c[0] = tmp[0];
    c[1..].copy_from_slice(&tmp[2..n + 1]);
    return 0;
}

#[no_mangle]
pub unsafe extern "C" fn make_cfft_plan(length: usize) -> cfft_plan {
    if length == 0 {
        return null_mut();
    }
    let plan: cfft_plan = Box::into_raw(Box::new(cfft_plan_i {
        packplan: null_mut(),
        blueplan: null_mut(),
    }));
    if plan.is_null() {
        return null_mut();
    }
    (*plan).blueplan = null_mut() as fftblue_plan;
    (*plan).packplan = null_mut() as cfftp_plan;
    if (length < 50) || (largest_prime_factor(length) <= ((length as f64).sqrt() as usize)) {
        (*plan).packplan = make_cfftp_plan(length);
        if (*plan).packplan.is_null() {
            let _ = Box::from_raw(plan);
            return null_mut();
        }
        return plan;
    }
    let comp1: f64 = cost_guess(length);
    let mut comp2: f64 = 2.0 * cost_guess(good_size(2 * length - 1));
    comp2 *= 1.5;
    if comp2 < comp1 {
        (*plan).blueplan = make_fftblue_plan(length);
        if (*plan).blueplan.is_null() {
            let _ = Box::from_raw(plan);
            return null_mut();
        }
    } else {
        (*plan).packplan = make_cfftp_plan(length);
        if (*plan).packplan.is_null() {
            let _ = Box::from_raw(plan);
            return null_mut();
        }
    }
    return plan;
}
#[no_mangle]
pub unsafe extern "C" fn destroy_cfft_plan(plan: cfft_plan) {
    if !(*plan).blueplan.is_null() {
        destroy_fftblue_plan((*plan).blueplan);
    }
    if !(*plan).packplan.is_null() {
        destroy_cfftp_plan((*plan).packplan);
    }
    let _ = Box::from_raw(plan);
}
#[no_mangle]
pub unsafe extern "C" fn cfft_backward(plan: cfft_plan, c: *mut f64, fct: f64) -> i32 {
    assert!(!plan.is_null(), "null");
    assert!(!c.is_null(), "null");
    if !(*plan).packplan.is_null() {
        let len = (*(*plan).packplan).length;
        let tmp_c = std::slice::from_raw_parts_mut(c, len);
        return cfftp_backward(&mut *(*plan).packplan, tmp_c, fct);
    }
    let n: usize = (*(*plan).blueplan).n;
    let tmp_c = from_raw_parts_mut(c, n * 2);
    return cfftblue_backward(&mut *(*plan).blueplan, tmp_c, fct);
}
#[no_mangle]
pub unsafe extern "C" fn cfft_forward(plan: cfft_plan, c: *mut f64, fct: f64) -> i32 {
    assert!(!plan.is_null(), "null");
    assert!(!c.is_null(), "null");
    if !(*plan).packplan.is_null() {
        let len = (*(*plan).packplan).length;
        let c = std::slice::from_raw_parts_mut(c, len);
        return cfftp_forward(&mut *(*plan).packplan, c, fct);
    }
    let n: usize = (*(*plan).blueplan).n;
    let tmp_c = from_raw_parts_mut(c, n * 2);
    return cfftblue_forward(&mut *(*plan).blueplan, tmp_c, fct);
}

#[no_mangle]
pub unsafe extern "C" fn make_rfft_plan(length: usize) -> rfft_plan {
    if length == 0 {
        return null_mut();
    }
    let tmp_rfft_plan_i = rfft_plan_i {
        packplan: null_mut(),
        blueplan: null_mut(),
    };
    let plan: rfft_plan = Box::into_raw(Box::new(tmp_rfft_plan_i));
    if plan.is_null() {
        return null_mut();
    }
    (*plan).blueplan = null_mut();
    (*plan).packplan = null_mut();

    let length_sqrt = (length as f64).sqrt() as usize;
    if (length < 50) || (largest_prime_factor(length) <= length_sqrt) {
        (*plan).packplan = make_rfftp_plan(length);
        if (*plan).packplan.is_null() {
            let _ = Box::from_raw(plan);
            return null_mut();
        }
        return plan;
    }
    let comp1: f64 = 0.5 * cost_guess(length);
    let mut comp2: f64 = 2.0 * cost_guess(good_size(2 * length - 1));
    comp2 *= 1.5;
    if comp2 < comp1 {
        (*plan).blueplan = make_fftblue_plan(length);
        if (*plan).blueplan.is_null() {
            let _ = Box::from_raw(plan);
            return null_mut();
        }
    } else {
        (*plan).packplan = make_rfftp_plan(length);
        if (*plan).packplan.is_null() {
            let _ = Box::from_raw(plan);
            return null_mut();
        }
    }
    return plan;
}

#[no_mangle]
pub unsafe extern "C" fn destroy_rfft_plan(plan: rfft_plan) {
    assert!(!plan.is_null(), "null");
    if !(*plan).blueplan.is_null() {
        destroy_fftblue_plan((*plan).blueplan);
    }
    if !(*plan).packplan.is_null() {
        destroy_rfftp_plan((*plan).packplan);
    }
    let _ = Box::from_raw(plan);
}
#[no_mangle]
pub unsafe extern "C" fn rfft_length(plan: rfft_plan) -> usize {
    assert!(!plan.is_null(), "null");
    if !(*plan).packplan.is_null() {
        return (*(*plan).packplan).length;
    }
    return (*(*plan).blueplan).n;
}
#[no_mangle]
pub unsafe extern "C" fn cfft_length(plan: cfft_plan) -> usize {
    assert!(!plan.is_null(), "null");
    if !(*plan).packplan.is_null() {
        return (*(*plan).packplan).length;
    }
    return (*(*plan).blueplan).n;
}
#[no_mangle]
pub unsafe extern "C" fn rfft_backward(plan: rfft_plan, c: *mut f64, fct: f64) -> i32 {
    assert!(!plan.is_null(), "null");
    assert!(!c.is_null(), "null");
    if !(*plan).packplan.is_null() {
        let n: usize = (*(*plan).packplan).length;
        let tmp_c = from_raw_parts_mut(c, n);
        return rfftp_backward(&mut *(*plan).packplan, tmp_c, fct);
    } else {
        let n: usize = (*(*plan).blueplan).n;
        let tmp_c = from_raw_parts_mut(c, n);
        return rfftblue_backward(&mut *(*plan).blueplan, tmp_c, fct);
    }
}
#[no_mangle]
pub unsafe extern "C" fn rfft_forward(plan: rfft_plan, c: *mut f64, fct: f64) -> i32 {
    assert!(!plan.is_null(), "null");
    assert!(!c.is_null(), "null");
    if !(*plan).packplan.is_null() {
        let n: usize = (*(*plan).packplan).length;
        let tmp_c = from_raw_parts_mut(c, n);
        return rfftp_forward(&mut *(*plan).packplan, tmp_c, fct);
    } else {
        let n: usize = (*(*plan).blueplan).n;
        let tmp_c = from_raw_parts_mut(c, n);
        return rfftblue_forward(&mut *(*plan).blueplan, tmp_c, fct);
    }
}

//=========================================================
//=========================================================
pub fn make_rfft_plan_safe(length: usize) -> RfftPlan {
    let mut plan: RfftPlan = RfftPlan {
        packplan: RfftpPlan {
            length: 0,
            nfct: 0,
            mem: Vec::new(),
            fct: [RfftpFctdata {
                fct: 0,
                tw_index: 0,
                tws_index: 0,
            }; NFCT],
        },
        blueplan: FftbluePlan {
            n: length,
            n2: good_size(length * 2 - 1),
            mem: Vec::new(),
            plan: CfftpPlan {
                length: 0,
                nfct: 0,
                mem: Vec::new(),
                fct: [CfftpFctdata {
                    fct: 0,
                    tw_index: 0,
                    tws_index: 0,
                }; NFCT],
            },
        },
    };
    //plan.blueplan.mem = vec![0.0f64; 2 * (plan.blueplan).n + 2 * (plan.blueplan).n2];
    let length_sqrt = (length as f64).sqrt() as usize;
    if (length < 50) || (largest_prime_factor(length) <= length_sqrt) {
        plan.packplan = make_rfftp_plan_safe(length);
        return plan;
    }
    let comp1: f64 = 0.5 * cost_guess(length);
    let mut comp2: f64 = 2.0 * cost_guess(good_size(2 * length - 1));
    comp2 *= 1.5;
    if comp2 < comp1 {
        plan.blueplan = make_fftblue_plan_safe(length);
    } else {
        plan.packplan = make_rfftp_plan_safe(length);
    }
    return plan;
}

fn make_fftblue_plan_safe(length: usize) -> FftbluePlan {
    let tmp_fct = [CfftpFctdata {
        fct: 0,
        tw_index: 0,
        tws_index: 0,
    }; NFCT];
    let tmp_cfftp_plan = CfftpPlan {
        length: 0,
        nfct: 0,
        mem: Vec::new(),
        fct: tmp_fct,
    };
    let mut plan: FftbluePlan = FftbluePlan {
        n: 0,
        n2: 0,
        plan: tmp_cfftp_plan,
        mem: Vec::new(),
    };

    (plan).n = length;
    let n2 = good_size(length * 2 - 1);
    (plan).n2 = n2;

    (plan).mem = vec![0.0f64; 2 * (plan).n + 2 * (plan).n2];

    let mut tmp = vec![0.0f64; 4 * length];
    sincos_2pibyn(2 * length, tmp.as_mut_slice());
    let bk = (plan).mem.as_mut_slice();
    bk[0] = 1.0;
    bk[1] = 0.0;
    let mut coeff: usize = 0;
    let mut m: usize = 1;
    while m < length {
        coeff = coeff + m * 2 - 1;
        if coeff >= length * 2 {
            coeff -= length * 2;
        }
        bk[m * 2] = tmp[2 * coeff];
        bk[m * 2 + 1] = tmp[2 * coeff + 1];
        m += 1;
    }
    let xn2: f64 = 1.0f64 / (plan).n2 as f64;
    let bkf_index = length * 2;
    bk[bkf_index] = bk[0] * xn2;
    bk[bkf_index + 1] = bk[1] * xn2;
    m = 2;
    while m < length * 2 {
        bk[bkf_index + (2 * n2 - m)] = bk[m] * xn2;
        bk[bkf_index + m] = bk[bkf_index + (2 * n2 - m)];

        bk[bkf_index + (2 * n2 - m + 1)] = bk[m + 1] * xn2;
        bk[bkf_index + m + 1] = bk[bkf_index + (2 * n2 - m + 1)];
        m += 2;
    }

    m = length * 2;
    while m <= (2 * ((plan).n2) - 2 * ((plan).n) + 1) {
        bk[bkf_index + m] = 0.0f64;
        m += 1;
    }
    (plan).plan = make_cfftp_plan_safe((plan).n2);
    if cfftp_forward_safe(&mut (plan).plan, &mut bk[bkf_index..], 1.0f64) != 0_i32 {
        return Default::default();
    }
    return plan;
}

fn cfftp_forward_safe(plan: &mut CfftpPlan, c: &mut [f64], fct: f64) -> i32 {
    let c = unsafe { std::mem::transmute::<&mut [f64], &mut [cmplx]>(c) };

    return pass_all_safe(plan, c, fct, -1);
}

fn cfftp_backward_safe(plan: &mut CfftpPlan, c: &mut [f64], fct: f64) -> i32 {
    let c = unsafe { std::mem::transmute::<&mut [f64], &mut [cmplx]>(c) };

    return pass_all_safe(plan, c, fct, 1);
}

fn pass_all_safe(plan: &mut CfftpPlan, c: &mut [cmplx], fct: f64, sign: i32) -> i32 {
    let len = plan.length;
    if len == 1 {
        return 0;
    }
    let mut l1 = 1;
    let nf = plan.nfct;
    let mut ch: Vec<cmplx> = vec![cmplx { r: 0.0, i: 0.0 }; len];

    let c_ptr = c.as_ptr();
    let mut p1: &mut [cmplx] = &mut c[..len];
    let mut p2: &mut [cmplx] = ch.as_mut_slice();
    let mut k1 = 0;
    while k1 < nf {
        let ip: usize = plan.fct[k1].fct;
        let l2: usize = ip * l1;
        let ido: usize = len / l2;
        let wa = &plan.mem[plan.fct[k1].tw_index..];
        if ip == 4 {
            if sign > 0 {
                pass4b(ido, l1, p1, p2, wa);
            } else {
                pass4f(ido, l1, p1, p2, wa);
            };
        } else if ip == 2 {
            if sign > 0 {
                pass2b(ido, l1, p1, p2, wa);
            } else {
                pass2f(ido, l1, p1, p2, wa);
            };
        } else if ip == 3 {
            if sign > 0 {
                pass3b(ido, l1, p1, p2, wa);
            } else {
                pass3f(ido, l1, p1, p2, wa);
            };
        } else if ip == 5 {
            if sign > 0 {
                pass5b(ido, l1, p1, p2, wa);
            } else {
                pass5f(ido, l1, p1, p2, wa);
            };
        } else if ip == 7 {
            pass7(ido, l1, p1, p2, wa, sign as i64);
        } else if ip == 11 {
            pass11(ido, l1, p1, p2, wa, sign as i64);
        } else {
            let csarr = &plan.mem[plan.fct[k1].tws_index..];
            let res = passg(ido, ip, l1, p1, p2, wa, csarr, sign as i64);
            if res != 0 {
                return -1;
            }
            std::mem::swap(&mut p1, &mut p2);
        }
        std::mem::swap(&mut p1, &mut p2);
        l1 = l2;
        k1 += 1;
    }
    if p1.as_ptr() != c_ptr {
        if fct != 1.0f64 {
            for i in 0..len {
                p2[i].r = p1[i].r * fct;
                p2[i].i = p1[i].i * fct;
            }
        } else {
            p2.copy_from_slice(p1);
        }
    } else if fct != 1.0f64 {
        for i in 0..len {
            p1[i].r *= fct;
            p1[i].i *= fct;
        }
    }
    return 0;
}

fn make_rfftp_plan_safe(length: usize) -> RfftpPlan {
    let mut plan: RfftpPlan = RfftpPlan {
        length,
        nfct: 0,
        mem: Vec::new(),
        fct: [RfftpFctdata {
            fct: 0,
            tw_index: 0,
            tws_index: 0,
        }; NFCT],
    };

    if length == 1 {
        return plan;
    }
    if rfftp_factorize_safe(&mut plan) != 0_i32 {
        return plan;
    }
    let tws = rfftp_twsize_safe(&mut plan);
    plan.mem = vec![0.0f64; tws];
    if rfftp_comp_twiddle_safe(&mut plan) != 0_i32 {
        return plan;
    }
    return plan;
}

fn rfftp_factorize_safe(plan: &mut RfftpPlan) -> i32 {
    let mut length = plan.length;
    let mut nfct: usize = 0;
    while (length % 4) == 0 {
        if nfct >= NFCT {
            return -1_i32;
        }
        let fresh5 = nfct;
        nfct += 1;
        plan.fct[fresh5].fct = 4;
        length >>= 2_i32
    }
    if (length % 2) == 0 {
        length >>= 1_i32;
        if nfct >= NFCT {
            return -1_i32;
        }
        let fresh6 = nfct;
        nfct += 1;
        plan.fct[fresh6].fct = 2;
        let tmp: usize = plan.fct[0].fct;
        plan.fct[0].fct = plan.fct[nfct - 1].fct;
        plan.fct[nfct - 1].fct = tmp;
    }
    let mut maxl: usize = ((length as f64).sqrt() as usize) + 1;
    let mut divisor: usize = 3;
    while length > 1 && divisor < maxl {
        if (length % divisor) == 0 {
            while (length % divisor) == 0 {
                if nfct >= NFCT {
                    return -1_i32;
                }
                let fresh7 = nfct;
                nfct += 1;
                plan.fct[fresh7].fct = divisor;
                length /= divisor;
            }
            maxl = ((length as f64).sqrt() as usize) + 1;
        }
        divisor += 2;
    }
    if length > 1 {
        let fresh8 = nfct;
        nfct += 1;
        plan.fct[fresh8].fct = length;
    }
    plan.nfct = nfct;
    return 0_i32;
}

fn rfftp_twsize_safe(plan: &mut RfftpPlan) -> usize {
    let mut twsize: usize = 0;
    let mut l1: usize = 1;
    for k in 0..plan.nfct {
        let ip = plan.fct[k].fct;
        let ido = plan.length / (l1 * ip);
        twsize += (ip - 1) * (ido - 1);
        if ip > 5 {
            twsize += 2 * ip;
        }
        l1 *= ip;
    }
    return twsize;
}

fn rfftp_comp_twiddle_safe(plan: &mut RfftpPlan) -> i32 {
    let length: usize = plan.length;
    let mut twid = vec![0.0f64; 2 * length];
    sincos_2pibyn_half(length, twid.as_mut_slice());
    let mut l1: usize = 1;
    let mut k: usize = 0;
    let nfct = plan.nfct;
    let mut index = 0;
    while k < nfct {
        let ip = plan.fct[k].fct;
        let ido: usize = length / (l1 * ip);
        if k < nfct - 1 {
            plan.fct[k].tw_index = index;
            let mut j: usize = 1;
            while j < ip {
                let mut i: usize = 1;
                let tw = &mut plan.mem[plan.fct[k].tw_index..];
                while i <= ((ido - 1) / 2) {
                    tw[(j - 1) * (ido - 1) + 2 * i - 2] = twid[2 * j * l1 * i];
                    tw[(j - 1) * (ido - 1) + 2 * i - 1] = twid[2 * j * l1 * i + 1];
                    i += 1;
                }
                j += 1;
            }
            index += (ip - 1) * (ido - 1);
        }
        if ip > 5 {
            plan.fct[k].tws_index = index;
            let tws = &mut plan.mem[index..];
            tws[0] = 1.0f64;
            tws[1] = 0.0f64;
            let mut i: usize = 1;
            while i <= (ip >> 1) {
                tws[2 * i] = twid[2 * i * (length / ip)];
                tws[2 * i + 1] = twid[2 * i * (length / ip) + 1];
                tws[2 * (ip - i)] = twid[2 * i * (length / ip)];
                tws[2 * (ip - i) + 1] = -twid[2 * i * (length / ip) + 1];
                i += 1;
            }
            index += 2 * ip;
        }
        l1 *= ip;
        k += 1;
    }
    return 0_i32;
}

pub fn rfft_forward_safe(plan: &mut RfftPlan, c: &mut [f64], fct: f64) -> i32 {
    if plan.packplan.length != 0 {
        return rfftp_forward_safe(&mut plan.packplan, c, fct);
    } else {
        return rfftblue_forward_safe(&mut plan.blueplan, c, fct);
    }
}

fn rfftp_forward_safe(plan: &mut RfftpPlan, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = plan.length;
    if n == 1 {
        return 0;
    }
    let mut l1: usize = n;
    let nf: usize = plan.nfct;
    let mut ch = vec![0.0f64; n];
    let mut p1: &mut [f64] = &mut c[..n];
    let mut p2: &mut [f64] = ch.as_mut_slice();
    let mut k1: usize = 0;
    while k1 < nf {
        let k: usize = nf - k1 - 1;
        let ip: usize = plan.fct[k].fct;
        let ido: usize = n / l1;
        let wa = &plan.mem[plan.fct[k].tw_index..];
        l1 /= ip;
        if ip == 4 {
            radf4(ido, l1, p1, p2, wa);
        } else if ip == 2 {
            radf2(ido, l1, p1, p2, wa);
        } else if ip == 3 {
            radf3(ido, l1, p1, p2, wa);
        } else if ip == 5 {
            radf5(ido, l1, p1, p2, wa);
        } else {
            let csarr = &plan.mem[plan.fct[k].tws_index..];
            radfg(ido, ip, l1, p1, p2, wa, csarr);
            std::mem::swap(&mut p1, &mut p2);
        }
        std::mem::swap(&mut p1, &mut p2);
        k1 += 1;
    }
    copy_and_norm(p2, p1, n, fct);
    return 0;
}

fn rfftblue_forward_safe(plan: &mut FftbluePlan, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = plan.n;
    let mut tmp = vec![0.0f64; 2 * n];
    let mut m: usize = 0;
    while m < n {
        tmp[2 * m] = c[m];
        tmp[2 * m + 1] = 0.0f64;
        m += 1;
    }
    let res = fftblue_fft_safe(plan, tmp.as_mut_slice(), -1, fct);
    if res != 0 {
        return -1;
    }
    c[0] = tmp[0];
    c[1..n].copy_from_slice(&tmp[2..n + 1]);
    return 0;
}

fn fftblue_fft_safe(plan: &mut FftbluePlan, c: &mut [f64], isign: i32, fct: f64) -> i32 {
    let n: usize = plan.n;
    let n2: usize = plan.n2;
    let mut akf = vec![0.0f64; 2 * n2];
    let bk = plan.mem.as_slice();
    let bkf_index = n * 2;
    if isign > 0 {
        let mut m: usize = 0;
        while m < (n * 2) {
            akf[m] = c[m] * bk[m] - c[m + 1] * bk[m + 1];
            akf[m + 1] = c[m] * bk[m + 1] + c[m + 1] * bk[m];
            m += 2;
        }
    } else {
        let mut m: usize = 0;
        while m < (n * 2) {
            akf[m] = c[m] * bk[m] + c[m + 1] * bk[m + 1];
            akf[m + 1] = -c[m] * bk[m + 1] + c[m + 1] * bk[m];
            m += 2;
        }
    }
    let mut m: usize = n * 2;
    while m < (n2 * 2) {
        akf[m] = 0.0;
        m += 1;
    }
    let res = cfftp_forward_safe(&mut plan.plan, akf.as_mut_slice(), fct);
    if res != 0 {
        return -1;
    }
    if isign > 0 {
        let mut m: usize = 0;
        while m < (n2 * 2) {
            let im: f64 = -akf[m] * bk[bkf_index + m + 1] + akf[m + 1] * bk[bkf_index + m];
            akf[m] = akf[m] * bk[bkf_index + m] + akf[m + 1] * bk[bkf_index + m + 1];
            akf[m + 1] = im;
            m += 2;
        }
    } else {
        let mut m: usize = 0;
        while m < (n2 * 2) {
            let im: f64 = akf[m] * bk[bkf_index + m + 1] + akf[m + 1] * bk[bkf_index + m];
            akf[m] = akf[m] * bk[bkf_index + m] - akf[m + 1] * bk[bkf_index + m + 1];
            akf[m + 1] = im;
            m += 2;
        }
    }
    let res = cfftp_backward_safe(&mut plan.plan, akf.as_mut_slice(), 1.0f64);
    if res != 0 {
        return -1;
    }
    if isign > 0 {
        let mut m: usize = 0;
        while m < (n * 2) {
            c[m] = bk[m] * akf[m] - bk[m + 1] * akf[m + 1];
            c[m + 1] = bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
            m += 2;
        }
    } else {
        let mut m: usize = 0;
        while m < (n * 2) {
            c[m] = bk[m] * akf[m] + bk[m + 1] * akf[m + 1];
            c[m + 1] = -bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
            m += 2;
        }
    }
    return 0;
}

fn make_cfftp_plan_safe(length: usize) -> CfftpPlan {
    let mut plan: CfftpPlan = Default::default();

    (plan).length = length;
    (plan).nfct = 0;

    if length == 1 {
        return plan;
    }
    if cfftp_factorize_safe(&mut plan) != 0_i32 {
        return Default::default();
    }
    let tws = cfftp_twsize_safe(&mut plan);
    (plan).mem = vec![cmplx { r: 0.0, i: 0.0 }; tws];

    if cfftp_comp_twiddle_safe(&mut plan) != 0_i32 {
        return Default::default();
    }
    return plan;
}

fn cfftp_factorize_safe(plan: &mut CfftpPlan) -> i32 {
    let mut length = plan.length;
    let mut nfct: usize = 0;
    while (length % 4) == 0 {
        if nfct >= NFCT {
            return -1;
        }
        let fresh1 = nfct;
        nfct += 1;
        plan.fct[fresh1].fct = 4;
        length >>= 2;
    }
    if (length % 2) == 0 {
        length >>= 1;
        if nfct >= NFCT {
            return -1;
        }
        let fresh2 = nfct;
        nfct += 1;
        plan.fct[fresh2].fct = 2;
        let tmp = plan.fct[0].fct;
        plan.fct[0].fct = plan.fct[nfct - 1].fct;
        plan.fct[nfct - 1].fct = tmp;
    }
    let mut maxl: usize = ((length as f64).sqrt() as usize) + 1;
    let mut divisor: usize = 3;
    while length > 1 && divisor < maxl {
        if (length % divisor) == 0 {
            while (length % divisor) == 0 {
                if nfct >= NFCT {
                    return -1;
                }
                let fresh3 = nfct;
                nfct += 1;
                plan.fct[fresh3].fct = divisor;
                length /= divisor;
            }
            maxl = ((length as f64).sqrt() as usize) + 1;
        }
        divisor += 2;
    }
    if length > 1 {
        let fresh4 = nfct;
        nfct += 1;
        plan.fct[fresh4].fct = length;
    }
    plan.nfct = nfct;
    return 0;
}

fn cfftp_twsize_safe(plan: &mut CfftpPlan) -> usize {
    let mut twsize: usize = 0;
    let mut l1: usize = 1;
    let mut k: usize = 0;
    let nfct = plan.nfct;
    while k < nfct {
        let ip: usize = plan.fct[k].fct;
        let ido: usize = plan.length / (l1 * ip);
        twsize += (ip - 1) * (ido - 1);
        if ip > 11 {
            twsize += ip;
        }
        l1 *= ip;
        k += 1;
    }
    return twsize;
}

fn cfftp_comp_twiddle_safe(plan: &mut CfftpPlan) -> i32 {
    let length: usize = plan.length;
    let mut twid = vec![0.0f64; 2 * length];
    sincos_2pibyn(length, twid.as_mut_slice());
    let mut l1: usize = 1;
    let mut memofs: usize = 0;
    let mut k: usize = 0;
    let nfct = plan.nfct;
    while k < nfct {
        let ip: usize = plan.fct[k].fct;
        let ido: usize = length / (l1 * ip);
        plan.fct[k].tw_index = memofs;
        memofs += (ip - 1) * (ido - 1);

        for j in 1..ip {
            let mut i: usize = 1;
            let tw = &mut plan.mem[plan.fct[k].tw_index..];
            while i < ido {
                tw[(j - 1) * (ido - 1) + i - 1].r = twid[2 * j * l1 * i];
                tw[(j - 1) * (ido - 1) + i - 1].i = twid[2 * j * l1 * i + 1];
                i += 1;
            }
        }
        if ip > 11 {
            plan.fct[k].tws_index = memofs;
            let tws = &mut plan.mem[plan.fct[k].tws_index..];
            for j in 0..ip {
                tws[j].r = twid[2 * j * l1 * ido];
                tws[j].i = twid[2 * j * l1 * ido + 1];
            }
            memofs += ip;
        }
        l1 *= ip;
        k += 1;
    }
    return 0_i32;
}

pub fn rfft_backward_safe(plan: &mut RfftPlan, c: &mut [f64], fct: f64) -> i32 {
    if plan.packplan.length != 0 {
        return rfftp_backward_safe(&mut plan.packplan, c, fct);
    } else {
        return rfftblue_backward_safe(&mut plan.blueplan, c, fct);
    }
}

fn rfftp_backward_safe(plan: &mut RfftpPlan, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = plan.length;
    if n == 1 {
        return 0;
    }
    let mut l1: usize = 1;
    let nf: usize = plan.nfct;
    let mut ch = vec![0.0f64; n];
    let c_ptr = c.as_mut_ptr();
    let mut p1: &mut [f64] = &mut c[..n];
    let mut p2: &mut [f64] = ch.as_mut_slice();
    let mut k: usize = 0;
    while k < nf {
        let ip: usize = plan.fct[k].fct;
        let ido: usize = n / (ip * l1);
        let tw = plan.mem.as_slice();
        let wa = &tw[plan.fct[k].tw_index..];
        if ip == 4 {
            radb4(ido, l1, p1, p2, wa);
        } else if ip == 2 {
            radb2(ido, l1, p1, p2, wa);
        } else if ip == 3 {
            radb3(ido, l1, p1, p2, wa);
        } else if ip == 5 {
            radb5(ido, l1, p1, p2, wa);
        } else {
            let tws = plan.mem.as_slice();
            let csarr = &tws[plan.fct[k].tws_index..];
            radbg(ido, ip, l1, p1, p2, wa, csarr);
        }
        std::mem::swap(&mut p1, &mut p2);
        l1 *= ip;
        k += 1;
    }
    unsafe {
        copy_and_norm(from_raw_parts_mut(c_ptr, n), p1, n, fct);
    }
    return 0;
}

fn rfftblue_backward_safe(plan: &mut FftbluePlan, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = plan.n;
    let mut tmp = vec![0.0f64; 2 * n];
    tmp[0] = c[0];
    tmp[1] = 0.0f64;
    tmp[2..n + 1].copy_from_slice(&c[1..n]);
    if n & 1 == 0 {
        tmp[n + 1] = 0.0f64;
    }
    let mut m: usize = 2;
    while m < n {
        tmp[2 * n - m] = tmp[m];
        tmp[2 * n - m + 1] = -tmp[m + 1];
        m += 2;
    }
    let res = fftblue_fft_safe(plan, tmp.as_mut_slice(), 1, fct);
    if res != 0 {
        return -1;
    }

    for m in 0..n {
        c[m] = tmp[2 * m];
    }
    return 0;
}

pub fn make_cfft_plan_safe(length: usize) -> CfftPlan {
    let tmp_packplan: CfftpPlan = Default::default();
    let tmp_blueplan: FftbluePlan = Default::default();
    let mut plan = CfftPlan {
        packplan: tmp_packplan,
        blueplan: tmp_blueplan,
    };
    if (length < 50) || (largest_prime_factor(length) <= ((length as f64).sqrt() as usize)) {
        (plan).packplan = make_cfftp_plan_safe(length);
        return plan;
    }
    let comp1: f64 = cost_guess(length);
    let mut comp2: f64 = 2.0 * cost_guess(good_size(2 * length - 1));
    comp2 *= 1.5;
    if comp2 < comp1 {
        (plan).blueplan = make_fftblue_plan_safe(length);
    } else {
        (plan).packplan = make_cfftp_plan_safe(length);
    }
    return plan;
}

pub fn cfft_forward_safe(plan: &mut CfftPlan, c: &mut [f64], fct: f64) -> i32 {
    if plan.packplan.length != 0 {
        return cfftp_forward_safe(&mut plan.packplan, c, fct);
    }
    return cfftblue_forward_safe(&mut plan.blueplan, c, fct);
}

fn cfftblue_forward_safe(plan: &mut FftbluePlan, c: &mut [f64], fct: f64) -> i32 {
    return fftblue_fft_safe(plan, c, -1, fct);
}

pub fn cfft_backward_safe(plan: &mut CfftPlan, c: &mut [f64], fct: f64) -> i32 {
    if plan.packplan.length != 0 {
        return cfftp_backward_safe(&mut plan.packplan, c, fct);
    }
    return cfftblue_backward_safe(&mut plan.blueplan, c, fct);
}

fn cfftblue_backward_safe(plan: &mut FftbluePlan, c: &mut [f64], fct: f64) -> i32 {
    return fftblue_fft_safe(plan, c, 1, fct);
}
