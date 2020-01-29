use std::fs::File;
use std::path::Path;
use std::path::PathBuf;
use std::env;
use std::error::Error;


#[test]
fn solve() -> Result<(), Box<dyn Error>> {
    use sparse21::System;

    let mut p = env::current_dir()?;
    p.push("tests/data/mat0.mat");

    let s = System::from_file(&p)?;
    let (mut mat, rhs) = s.split();

    println!("Solving");
    let res = mat.solve(rhs)?;
    println!("Result: {:?}", res);

    return Ok(());
}

fn main() -> Result<(), Box<dyn Error>> {
    return solve();
}

