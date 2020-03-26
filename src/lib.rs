#[macro_use]
extern crate libc;
extern crate futures;
extern crate lapack;

use crate::futures::StreamExt;
use core::arch::x86_64::{
    __m128d, _mm_cvtsd_f64, _mm_mul_pd, _mm_set1_pd, _mm_set_pd, _mm_shuffle_pd,
};
use futures::{join, stream};
use libc::c_double;
use std::convert::TryInto;
use std::mem;
use std::{slice, vec};

// Things to do
// rename things, check that c header lib is also compat
// extract unsafe stuff
// Fill in uprh matrix and pass it to cholesky
// Compile accelerated and non using github.
// Send to Zhielang
//
struct LocalStiffness {
    matrix: Vec<__m128d>,
    index: usize,
}

#[no_mangle]
pub extern "C" fn solve(
    element_count: u16,
    element_nodes: u16,
    integration_points: u16,
    node_count: u16,
    plane_stress: u16,
    conditions_count: u16,
    reactions_count: u16,
    nu: c_double,
    E: c_double,
    body_force: c_double,
    nodes: *const c_double,
    elements: *const u16,
    reactions: *const u16,
    displacement_nodes: *const u16,
    displacement_given: *const c_double,
    displacement_result: *mut c_double,
    reaction_result: *mut c_double,
) {
    let size: usize = elements as usize;
    let (nodes, elements, reactions, displacement_nodes, displacement_given) = unsafe {
        let nodes = (1..node_count)
            .map(|i| {
                _mm_set_pd(
                    *nodes.offset((i * 2 + 1) as isize),
                    *nodes.offset((i * 2) as isize),
                )
            })
            .collect::<Vec<__m128d>>();

        let elements = (1..element_count)
            .map(|i| {
                slice::from_raw_parts(
                    elements.offset((i * element_nodes) as isize),
                    element_nodes as usize,
                )
                .to_vec()
            })
            .collect::<Vec<Vec<u16>>>();

        let reactions = (1..reactions_count)
            .map(|i| *reactions.offset(i as isize))
            .collect::<Vec<u16>>();

        let displacement_nodes = (1..conditions_count)
            .map(|i| slice::from_raw_parts(displacement_nodes.offset((i * 2) as isize), 2).to_vec())
            .collect::<Vec<Vec<u16>>>();

        let displacement_given = (1..conditions_count)
            .map(|i| *displacement_given.offset(i as isize))
            .collect::<Vec<f64>>();

        (
            nodes,
            elements,
            reactions,
            displacement_nodes,
            displacement_given,
        )
    };

    let equations = (nodes.len() * (nodes.len() + 1)) - displacement_nodes.len();
    // Construct Stiffness Matrix UPRH
    // map -> vec[f64] : c_uint
    // collect vec[f64] -> vec[vec[f64]]
    // fold vec[vec[f64]] -> vec[f64]
    // then move into memory
    let mut stiffness = (0..equations).map(|i| i as f64).collect::<Vec<f64>>();
    let stiffness_stream = stream::iter(elements)
        .map(|v| {
            // For each node, we produce a local stiffness matrix
            v.iter()
                .map(|index| {
                    let i: usize = (*index).try_into().unwrap();
                    let matrix = vec![nodes[i], _mm_mul_pd(nodes[i], _mm_set1_pd(2.0))];
                    LocalStiffness {
                        matrix,
                        index: *index as usize,
                    }
                })
                .collect::<Vec<LocalStiffness>>();
        })
        .collect::<Vec<_>>()
        .select_all()
        .map(|ke| {
            stiffness[ke.index] = _mm_cvtsd_f64(ke.matrix[0]);
            let temp = _mm_shuffle_pd(ke.matrix[0], ke.matrix[0], 1);
            stiffness[ke.index * 2 + 1] = _mm_cvtsd_f64(temp);
        });

    // Construct force vector and map to the provided solution vector.
    // then move into memory
    // let force_stream = stream::iter_ok::<_, ()>(nodes).map().fold().map();
    let forces = vec![1.0; equations];

    async {
        join!(stiffness_stream);
        // !join(stiffness_stream, force_stream);
    };

    // Do cholesky decomposition mkl::potrs
    unsafe {
        *reaction_result = mem::transmute_copy(&forces);
        lapack::c::dpotrf(
            lapack::c::Layout::RowMajor, //LAPACK_ROW_MAJOR
            b'U',
            equations as i32,
            &mut stiffness.as_slice(),
            1,
        );
        lapack::c::dpotrs(
            lapack::c::Layout::RowMajor, //LAPACK_ROW_MAJOR
            b'U',
            equations as i32,
            1,
            &stiffness.as_slice(),
            1,
            &mut forces.as_slice(),
            1,
        );
        *displacement_result = mem::transmute_copy(&forces);
    }

    // Solve the system with the Cholesky demcoposition mkl::potrs

    // unsafe {
    //    *c_double.offset(i as isize) = c[i];
    // }
}

fn safe_solve() {}

#[test]
fn it_works() {
    let a: Vec<f64> = vec![1f64, 2f64, 3f64];
    let b: Vec<f64> = vec![3f64, 2f64, 1f64];
    let c: Vec<f64> = multiply_safe(a, b);
    let expected: Vec<f64> = vec![3f64, 4f64, 3f64];

    assert!(c.len() == expected.len());

    for i in 0..c.len() {
        assert!(c[i] == expected[i]);
    }
}
