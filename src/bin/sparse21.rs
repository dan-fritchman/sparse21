use std::fs::File;
use std::path::Path;
use std::path::PathBuf;
use std::env;
use std::error::Error;


fn solve(p: &Path) -> Result<(), Box<dyn Error>> {
    use sparse21::System;
    let s = System::from_file(p)?;

    println!("Solving");
    let res = s.solve();
    println!("Result: {:?}", res);

    return Ok(());
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let p = Path::new(&args[1]);
    return solve(p);
}
