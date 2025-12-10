// change bias model easily here 
use ndarray::Array1;
use rustc_hash::FxHashMap;
use hdf5::File as hdf5_file;

pub fn read_bias_factors(h5_path: &str) -> FxHashMap<String, Vec<f32>> {
    let file = hdf5_file::open(h5_path).expect("Unable to open HDF5 file");
    let mut bias_factors: FxHashMap<String, Vec<f32>>= FxHashMap::default();
    // define range of chrom to use // throw the rest to save memory 
    let chromosomes: Vec<String> = (1..=22).map(|i| format!("chr{}", i)).collect();
    for key in chromosomes {
        let dataset = file.dataset(&key).expect("Unable to open dataset");
        // Read the dataset into an ndarray Array
        let array: Array1<f32> = dataset.read().expect("Unable to read dataset");
        // Convert the ndarray Array into a Vec<f32>
        let bias: Vec<f32> = array.iter().copied().collect();
        // let bias: Vec<f32> = dataset.read().expect("Unable to read dataset");
        bias_factors.insert(key, bias);   
    }
    bias_factors
}


pub fn cal_sum_bias(chrom: &String, start: usize, end: usize, bias_factors: &FxHashMap<String, Vec<f32>>) -> f32 {
    if let Some(bias) = bias_factors.get(chrom) {
        bias[start..end].iter().copied().sum() // NON END INCLUSIVE
    } else {
        0.167 * ((end - start) as f32)
    }
}

pub fn cal_mean_bias(chrom: &String, start: usize, end: usize, bias_factors: &FxHashMap<String, Vec<f32>>) -> f32 {
    if let Some(bias) = bias_factors.get(chrom) {
        let bias_sum: f32 = bias[start..end].iter().copied().sum();
        let bias_mean = bias_sum/((end-start) as f32);
        // 0.167_f32.max(bias_mean)
        bias_mean
    } else {
        //0.167
        0.0
    }
}
pub fn cal_max_bias(chrom: &String, start: usize, end: usize, bias_factors: &FxHashMap<String, Vec<f32>>) -> f32 {
    if let Some(bias) = bias_factors.get(chrom) {
        let bias_max = bias[start..=end].iter().copied().fold(0.0, f32::max); 
        0.167_f32.max(bias_max)
    } else {
        0.167 
    }
}

// Helper Functions 
pub fn calculate_mean(data: &Vec<f32>) -> f32 {
    data.iter().sum::<f32>()  / data.len() as f32
}

pub fn calculate_std(data: &Vec<f32>) -> f32 {
    let data_mean = calculate_mean(data);
    let variance = data.iter().map(|value| {
        let diff = data_mean - (*value as f32);
        diff * diff
    }).sum::<f32>() / data.len() as f32;

    variance.sqrt()
}
