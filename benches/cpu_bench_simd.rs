// use binius_rust::simd::binary_ntt_cache_gfni::WI_EVAL_CACHE;
// use binius_rust::simd::pcs::{commit, prove, verifier};
// #[cfg(target_arch = "x86_64")]
// use core::arch::x86_64::*;
// use criterion::{black_box, criterion_group, criterion_main, Criterion};

// #[cfg(all(target_arch = "x86_64", target_feature = "gfni"))]
// unsafe fn test_gfni(a: __m128i, b: __m128i) -> __m128i {
//     _mm_gf2p8mul_epi8(a, b)
// }
// fn detect_gfni() -> bool {
//     println!("Detecting GFNI support...");
//     #[cfg(target_arch = "x86_64")]
//     {
//         if is_x86_feature_detected!("gfni") {
//             println!("GFNI detected and enabled.");
//             true
//         } else {
//             println!("GFNI not supported.");
//             false
//         }
//     }
//     #[cfg(not(target_arch = "x86_64"))]
//     {
//         println!("Not running on x86_64 architecture.");
//         false
//     }
// }

// fn benchmark_commit(c: &mut Criterion) {
//     let evaluations = vec![1; 1 << 20];
//     // let evaluations = vec![1; 1 << 24];

//     {
//         let _unused = WI_EVAL_CACHE.lock().unwrap();
//     }

//     c.bench_function("commit_benchmark", |b| {
//         b.iter(|| {
//             let result = commit(black_box(&evaluations));
//             black_box(result);
//         })
//     });
// }

// fn benchmark_prove(c: &mut Criterion) {
//     let evaluations = vec![1; 1 << 20];
//     let evaluation_point = vec![1; 23];
//     // let evaluations = vec![1; 1 << 24];
//     // let evaluation_point = vec![1; 27];

//     {
//         let _unused = WI_EVAL_CACHE.lock().unwrap();
//     }
//     let commitment = commit(&evaluations);

//     c.bench_function("prove_benchmark", |b| {
//         b.iter(|| {
//             let result = prove(
//                 black_box(&commitment),
//                 black_box(&evaluations),
//                 black_box(&evaluation_point),
//             );
//             black_box(result);
//         })
//     });
// }

// fn benchmark_verifier(c: &mut Criterion) {
//     let evaluations = vec![1; 1 << 20];
//     let evaluation_point = vec![1; 23];
//     // let evaluations = vec![1; 1 << 24];
//     // let evaluation_point = vec![1; 27];

//     {
//         let _unused = WI_EVAL_CACHE.lock().unwrap();
//     }
//     let commitment = commit(&evaluations);
//     let proof = prove(&commitment, &evaluations, &evaluation_point);

//     c.bench_function("verifier_benchmark", |b| {
//         b.iter(|| {
//             let result = verifier(
//                 black_box(&commitment),
//                 black_box(&proof),
//                 black_box(&evaluation_point),
//             );
//             black_box(result);
//         })
//     });
// }

// criterion_group!(
//     benches,
//     benchmark_commit,
//     benchmark_prove,
//     benchmark_verifier
// );
// criterion_main!(benches);

use binius_rust::simd::binary_ntt_cache_gfni::WI_EVAL_CACHE;
use binius_rust::simd::pcs::{commit, prove, verifier};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
#[cfg(target_arch = "x86_64")]
use std::arch::is_x86_feature_detected;

#[cfg(target_arch = "x86_64")]
fn print_simd_features() {
    println!(
        "SIMD Features: GFNI: {}, SSE2: {}, AVX2: {}",
        is_x86_feature_detected!("gfni"),
        is_x86_feature_detected!("sse2"),
        is_x86_feature_detected!("avx2")
    );
}

fn benchmark_group_1(c: &mut Criterion) {
    let evaluations = vec![1; 1 << 20];
    let evaluation_point = vec![1; 23];

    {
        let _unused = WI_EVAL_CACHE.lock().unwrap();
    }

    let commitment = commit(&evaluations);

    c.bench_function("group_1_commit", |b| {
        b.iter(|| {
            let result = commit(black_box(&evaluations));
            black_box(result);
        })
    });

    c.bench_function("group_1_prove", |b| {
        b.iter(|| {
            let result = prove(
                black_box(&commitment),
                black_box(&evaluations),
                black_box(&evaluation_point),
            );
            black_box(result);
        })
    });

    c.bench_function("group_1_verifier", |b| {
        let proof = prove(&commitment, &evaluations, &evaluation_point);
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

fn benchmark_group_2(c: &mut Criterion) {
    let evaluations = vec![1; 1 << 24];
    let evaluation_point = vec![1; 27];

    {
        let _unused = WI_EVAL_CACHE.lock().unwrap();
    }

    let commitment = commit(&evaluations);

    c.bench_function("group_2_commit", |b| {
        b.iter(|| {
            let result = commit(black_box(&evaluations));
            black_box(result);
        })
    });

    c.bench_function("group_2_prove", |b| {
        b.iter(|| {
            let result = prove(
                black_box(&commitment),
                black_box(&evaluations),
                black_box(&evaluation_point),
            );
            black_box(result);
        })
    });

    c.bench_function("group_2_verifier", |b| {
        let proof = prove(&commitment, &evaluations, &evaluation_point);
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

fn benchmark_group_3(c: &mut Criterion) {
    let evaluations = vec![1; 1 << 28];
    let evaluation_point = vec![1; 31];

    {
        let _unused = WI_EVAL_CACHE.lock().unwrap();
    }

    let commitment = commit(&evaluations);

    c.bench_function("group_3_commit", |b| {
        b.iter(|| {
            let result = commit(black_box(&evaluations));
            black_box(result);
        })
    });

    c.bench_function("group_3_prove", |b| {
        b.iter(|| {
            let result = prove(
                black_box(&commitment),
                black_box(&evaluations),
                black_box(&evaluation_point),
            );
            black_box(result);
        })
    });

    c.bench_function("group_3_verifier", |b| {
        let proof = prove(&commitment, &evaluations, &evaluation_point);
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

fn benchmark_with_simd_detection(c: &mut Criterion) {
    #[cfg(target_arch = "x86_64")]
    print_simd_features(); // Print SIMD feature detection at the start
    benchmark_group_1(c);
    benchmark_group_2(c);
    benchmark_group_3(c);
}

criterion_group!(benches, benchmark_with_simd_detection);
criterion_main!(benches);
