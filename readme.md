
# sparse21

Solving large systems of linear equations using sparse matrix methods. 

[![Docs](https://docs.rs/sparse21/badge.svg)](docs.rs/sparse21)

```rust 
let mut m = sparse21::Matrix::from_entries(vec![
            (0, 0, 1.0),
            (0, 1, 1.0),
            (0, 2, 1.0),
            (1, 1, 2.0),
            (1, 2, 5.0),
            (2, 0, 2.0),
            (2, 1, 5.0),
            (2, 2, -1.0),
        ]);

        let soln = m.solve(vec![6.0, -4.0, 27.0]); 
        // => vec![5.0, 3.0, -2.0]
```

Sparse methods are primarily valuable for systems in which the number of non-zero entries is substantially less than the overall size of the matrix. Such situations are common in physical systems, including electronic circuit simulation. All elements of a sparse matrix are assumed to be zero-valued unless indicated otherwise. 

## Usage 

Sparse21 exposes two primary data structures: 

* `Matrix` represents an `f64`-valued sparse matrix
* `System` represents a system of linear equations of the form `Ax=b`, including a `Matrix` (A) and right-hand-side `Vec` (b).

Once matrices and systems have been created, their primary public method is `solve`, which returns a (dense) `Vec` solution-vector.


## Matrix 

Sparse21 matrices can be constructed from a handful of data-sources

`Matrix::new` creates an empty matrix, to which elements can be added via the `add_element` and `add_elements` methods. 

```rust
let mut m = Matrix::new();

m.add_element(0, 0, 11.0);
m.add_element(7, 0, 22.0);
m.add_element(0, 7, 33.0);
m.add_element(7, 7, 44.0);
```

```rust
let mut m = Matrix::new();

m.add_elements(vec![(0, 0, 11.0),
    (7, 0, 22.0),
    (0, 7, 33.0),
    (7, 7, 44.0)
]);
```

The arguments to `add_element` are a row (`usize`), column (`usize`), and value (`f64`). 
Adding elements (plural) via `add_elements` takes a vector of `(usize, usize, f64)` tuples, representing the row, col, and val. 

Unlike common mathematical notation, all locations in `sparse21` matrices and vectors are zero-indexed. 
Adding a non-zero at the "first" matrix element therefore implies calling `add_element(0, 0, val)`. 

Creating a `Matrix` from data entries with `Matrix::from_entries`: 

```rust
let mut m = Matrix::from_entries(vec![
            (2, 2, -1.0),
            (2, 1, 5.0),
            (2, 0, 2.0),
            (1, 2, 5.0),
            (1, 1, 2.0),
            (0, 2, 1.0),
            (0, 1, 1.0),
            (0, 0, 1.0),
        ]);
``` 

The `Matrix::identity` method returns a new identity matrix of size (n x n):

```rust
let mut m = Matrix::identity(3);
let soln = m.solve(vec![11.1, 30.3, 99.9]);
assert_eq!(soln, vec![11.1, 30.3, 99.9]);
```

### Solving 

Sparse21 matrices are built for solving equation-systems. The primary public method of a `Matrix` is `solve()`, which accepts a `Vec` right-hand-side as its sole argument, and returns a solution `Vec` of the same size. 
 
### Matrix Mutability 

You may have noticed all examples to date declare matrices as `mut`, perhaps unnecessarily. This is on purpose. The `Matrix::solve` method (un-rustily) modifies the matrix *in-place*. For larger matrices, the in-place modification saves orders of magnitude of memory, as well as time creating and destroying elements. While in-place self-modification falls out of line with the Rust ethos, it follows a long lineage of scientific computing tools for this and similar tasks. 

So: in order to be solved, matrices must be declared `mut`.   

## System 

Sparse21 equation-systems can be loaded directly from a (custom-format) file including the elements and RHS. 

```rust
let s = sparse21::System::from_file(Path::new("my_matrix.mat"))?;
let soln = s.solve();
```

It is often more convenient to split a `sparse21::System` into its constituent matrix and RHS. The `split` method returns a tuple of the two: 

```rust
let s = sparse21::System::from_file(&p)?;

// Split into parts, and solve 
// More or less equivalent to `s.solve()` 
let (mut mat, rhs) = s.split();
let res = mat.solve(rhs);
```

