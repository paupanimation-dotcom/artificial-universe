Artificial Universe — Phase 5i, step 1 delta
============================================

DIVISION, AND HEREDITY WITHOUT A GENOME.

5 files. 2 new, 3 modified. Extract over your project folder and overwrite.
Apply Phase 5h first if you have not.

  NEW
    crates/au-chem/src/vesicle.rs      geometry, fate, stochastic partition
    crates/au-chem/tests/vesicle.rs    15 tests

  MODIFIED
    CHANGELOG.md                       Phase 5i step 1 entry
    README.md                          status + Phase 5i section
    crates/au-chem/src/lib.rs          module registration + exports

VERIFIED
  cargo test --release  ->  230 passed, 0 failed  (was 215)
  Built from a clean tree with rustc 1.75; zero third-party dependencies.

NOTHING IS WIRED IN YET
  This step is pure au-chem. No protocell exists as an object, there is no
  storage and no lineage, and every world runs byte-identically to Phase 5h.
  What it establishes is WHEN a bag divides and WHAT its daughters get.
  Step 2 is entities: identity, contents, parentage, and the codec extension
  that lets a variable-size composition survive a save.

THE TWO ANCHORS
  v = 1/sqrt(2)   the reduced volume at which a membrane's area is exactly
                  enough to close two equal spheres. Derived, not chosen.
  1/(2*sqrt(n))   the relative error with which n copies are transmitted.
                  Fidelity is copy number; nothing sets it.
