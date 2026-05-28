use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bbduk(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bbduk");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    // bbduk golden may be in a different location; use the nearest fastq from preprocessing
    let fq = manifest
        .parent()
        .unwrap()
        .join("rsomics-fastp/tests/golden/test_r1.fastq");
    c.bench_function("rsomics-bbduk golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["--in1", fq.to_str().unwrap(), "--out1", "/dev/null"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bbduk);
criterion_main!(benches);
