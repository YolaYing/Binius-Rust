//! This benchmark is used to record the vanilla version of Binius's performance.

use binius_rust::binary_ntt_cache::WI_EVAL_CACHE;
use binius_rust::prover::prover;
use binius_rust::verifier::verifier;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn benchmark_prover(c: &mut Criterion) {
    // evaluations的大小2^20 u8, evaluations的数据类型是Vec<u8>
    // evaluation_point的数量是log2(len(evaluations) * 8) = log2(2^20 * 8) = 23
    let evaluations = vec![1; 1 << 20];
    let evaluation_point = vec![1; 23];

    {
        let _unused = WI_EVAL_CACHE.lock().unwrap();
    }

    c.bench_function("prover_benchmark", |b| {
        b.iter(|| {
            let result = prover(black_box(&evaluations), black_box(evaluation_point.clone()));
            black_box(result);
        })
    });
}

fn benchmark_verifier(c: &mut Criterion) {
    let evaluations = vec![1; 1 << 20];
    let evaluation_point = vec![1; 23];

    {
        let _unused = WI_EVAL_CACHE.lock().unwrap();
    }

    c.bench_function("verifier_benchmark", |b| {
        b.iter(|| {
            let proof = prover(black_box(&evaluations), black_box(evaluation_point.clone()));
            let result = verifier(black_box(proof));
            black_box(result);
        })
    });
}

criterion_group!(benches, benchmark_prover, benchmark_verifier);
criterion_main!(benches);
