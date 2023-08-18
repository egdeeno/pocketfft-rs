
use pocketfft_rs::*;


fn fill_random(data: &mut [f64], length: usize){
    use libc::rand;
    for m in 0..length{
        data[m] = (unsafe{rand()} as f64)/(i32::MAX as f64+1.0) - 0.5;
    }
}

fn errcalc(data: &[f64], odata: &[f64], length: usize )->f64{
    let mut sum:f64 = 0.0;
    let mut errsum:f64 = 0.0;
    for m in 0..length {
        errsum += (data[m]-odata[m])*(data[m]-odata[m]);
        sum += odata[m]*odata[m];
    }
    return (errsum/sum).sqrt();
}

const MAXLEN: usize = 8192;
const EPSILON: f64 = 2e-15;

#[test]
fn test_real_unsafe() {
    let mut data: [f64; MAXLEN];
    let mut odata: [f64; MAXLEN] = [0.0; MAXLEN];
    fill_random(&mut odata, MAXLEN);
    let mut errsum = 0.0;
    for length in 1..=MAXLEN {
        data = odata;
        unsafe {
            let plan = make_rfft_plan (length);
            rfft_forward (plan, data.as_mut_ptr(), 1.0);
            rfft_backward (plan, data.as_mut_ptr(), 1.0/length as f64);
            destroy_rfft_plan (plan);
        }
        let err = errcalc(&data, &odata, length);
        assert!(err<=EPSILON, "problem at real length {}: {}",length,err);
        errsum+=err;
    }
    println!("test_real errsum: {:e}",errsum);
}


#[test]
fn test_complex_unsafe() {
    let mut data:[f64; MAXLEN*2];
    let mut odata:[f64; MAXLEN*2] = [0.0; MAXLEN*2];
    fill_random(&mut odata, MAXLEN*2);
    let mut errsum = 0.0;
    for length in 1..=MAXLEN {
        data = odata;
        unsafe {
            let plan = make_cfft_plan (length);
            cfft_forward (plan, data.as_mut_ptr(), 1.0);
            cfft_backward (plan, data.as_mut_ptr(), 1.0/length as f64);
            destroy_cfft_plan (plan);
        }
        let err = errcalc(&data, &odata, length*2);
        assert!(err<=EPSILON, "problem at real length {}: {}",length,err);
        errsum+=err;
    }
    println!("test_complex errsum: {:e}",errsum);
}

#[test]
fn test_real_safe() {
    let mut data: [f64; MAXLEN];
    let mut odata: [f64; MAXLEN] = [0.0; MAXLEN];
    fill_random(&mut odata, MAXLEN);
    let mut errsum = 0.0;
    for length in 1..=MAXLEN {
        data = odata;
        let mut plan = make_rfft_plan_safe (length);
        rfft_forward_safe(&mut plan, data.as_mut_slice(), 1.0);
        rfft_backward_safe(&mut plan, data.as_mut_slice(), 1.0/length as f64);
        let err = errcalc(&data, &odata, length);
        assert!(err<=EPSILON, "problem at real length {}: {}",length,err);
        errsum+=err;
    }
    println!("test_real errsum: {:e}",errsum);
}

#[test]
fn test_complex_safe() {
    let mut data:[f64; MAXLEN*2];
    let mut odata:[f64; MAXLEN*2] = [0.0; MAXLEN*2];
    fill_random(&mut odata, MAXLEN*2);
    let mut errsum = 0.0;
    for length in 1..=MAXLEN {
        data = odata;
        let mut plan = make_cfft_plan_safe(length);
        cfft_forward_safe(&mut plan, data.as_mut_slice(), 1.0);
        cfft_backward_safe(&mut plan, data.as_mut_slice(), 1.0/length as f64);
        let err = errcalc(&data, &odata, length*2);
        assert!(err<=EPSILON, "problem at real length {}: {}",length,err);
        errsum+=err;
    }
    println!("test_complex errsum: {:e}",errsum);
}
