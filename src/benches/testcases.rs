use std::error::Error;


//#[test]
fn solve_all() -> Result<(), Box<dyn Error>> {
    use glob::glob;
    use sparse21::System;

    for f in glob("tests/data/*.mat")? {
        match f {
            Ok(f) => {
                let s = match System::from_file(&f) {
                    Err(e) => {
                        println!("Error: {}", e);
                        continue;
                    }
                    Ok(s) => s,
                };

                let (mut mat, rhs) = s.split();

//                println!("Factorizing");
//                mat.lu_factorize();
//                println!("Factorized");
//
//                if rhs.len() > 0 {
//                    println!("Solving");
//                    let res = mat.solve(rhs);
//                    println!("Result: {:?}", res);
//                }
            }
            Err(e) => println!("{:?}", e),
        }
    }
    return Ok(());
}
