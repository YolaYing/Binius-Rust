use binius_rust::binary_ntt_cache::WI_EVAL_CACHE;
use binius_rust::pcs::{commit, prove, verifier};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn benchmark_commit(c: &mut Criterion) {
    let evaluations = vec![1; 1 << 20];

    {
        let _unused = WI_EVAL_CACHE.lock().unwrap();
    }

    c.bench_function("commit_benchmark", |b| {
        b.iter(|| {
            let result = commit(black_box(&evaluations));
            black_box(result);
        })
    });
}

fn benchmark_prove(c: &mut Criterion) {
    let evaluations = vec![1; 1 << 20];
    let evaluation_point = vec![1; 23];

    {
        let _unused = WI_EVAL_CACHE.lock().unwrap();
    }
    let commitment = commit(&evaluations);

    c.bench_function("prove_benchmark", |b| {
        b.iter(|| {
            let result = prove(
                black_box(&commitment),
                black_box(&evaluations),
                black_box(&evaluation_point),
            );
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
    let commitment = commit(&evaluations);
    let proof = prove(&commitment, &evaluations, &evaluation_point);

    c.bench_function("verifier_benchmark", |b| {
        b.iter(|| {
            let result = verifier(
                black_box(&commitment),
                black_box(&proof),
                black_box(&evaluation_point),
            );
            black_box(result);
        })
    });
}

criterion_group!(
    benches,
    benchmark_commit,
    benchmark_prove,
    benchmark_verifier
);
criterion_main!(benches);
