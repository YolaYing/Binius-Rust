use criterion::{black_box, criterion_group, criterion_main, Criterion};

struct MyStruct {
    value: usize,
}

impl MyStruct {
    fn bit_length_method1(&self) -> usize {
        1 << (self.value.next_power_of_two().trailing_zeros() as usize)
    }

    fn bit_length_method2(&self) -> usize {
        let value_bit_length = self.value.bit_length();
        1 << (value_bit_length - 1).bit_length()
    }
}

fn benchmark_method1(c: &mut Criterion) {
    let my_struct = MyStruct { value: 10 };
    c.bench_function("bit_length_method1", |b| b.iter(|| my_struct.bit_length_method1()));
}

fn benchmark_method2(c: &mut Criterion) {
    let my_struct = MyStruct { value: 10 };
    c.bench_function("bit_length_method2", |b| b.iter(|| my_struct.bit_length_method2()));
}

criterion_group!(benches, benchmark_method1, benchmark_method2);
criterion_main!(benches);