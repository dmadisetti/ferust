#[macro_use]
extern crate libc;
extern crate futures;
extern crate itertools;
extern crate lapack_src;
extern crate lapacke;

use crate::futures::StreamExt;
use core::arch::x86_64::{
    __m128d, _mm_cvtsd_f64, _mm_mul_pd, _mm_set1_pd, _mm_set_pd, _mm_shuffle_pd,
};
use futures::{join, stream};
use itertools::Itertools;
use libc::c_double;
use std::collections::HashSet;
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
    // Cast from raw pointers into safer rust types.
    let element_nodes: usize = element_nodes as usize;
    let (nodes, elements, reactions, displacement_nodes, displacement_given): (
        std::vec::Vec<std::arch::x86_64::__m128d>,
        std::vec::Vec<std::vec::Vec<u16>>,
        std::vec::Vec<u16>,
        std::collections::HashSet<std::vec::Vec<u16>>,
        std::vec::Vec<f64>,
    ) = unsafe {
        let nodes = (0..node_count)
            .map(|i| {
                _mm_set_pd(
                    *nodes.offset((i * 2 + 1) as isize),
                    *nodes.offset((i * 2) as isize),
                )
            })
            .collect::<Vec<__m128d>>();

        let elements = (0..element_count)
            .map(|i| {
                slice::from_raw_parts(
                    elements.offset((i as isize) * element_nodes as isize),
                    element_nodes,
                )
                .to_vec()
            })
            .collect::<Vec<Vec<u16>>>();

        let reactions = (0..reactions_count)
            .map(|i| *reactions.offset(i as isize))
            .collect::<Vec<u16>>();

        let displacement_nodes = (0..conditions_count)
            .map(|i| slice::from_raw_parts(displacement_nodes.offset((i * 2) as isize), 2).to_vec())
            .collect::<HashSet<Vec<u16>>>();

        let displacement_given = (0..conditions_count)
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

    // TODO: Move everything into own module. Just here for proof of concept.

    let equations = nodes.len() * 2 - displacement_given.len();
    let stiffness_entries = equations * (equations + 1) / 2;
    // Construct Stiffness Matrix UPRH
    // map grid -> LocalStiffness
    // gather LocalStiffness -> Global Stiffness
    // then move into memory
    let mut stiffness = (0..stiffness_entries)
        .map(|i| i as f64)
        .collect::<Vec<f64>>();
    let stiffness_stream =
        stream::iter((0..(element_count - 1)).cartesian_product(0..(element_nodes - 1)))
            .map(|xy| {
                let index = elements[xy.0 as usize][xy.1 as usize];
                let i: usize = index.try_into().unwrap();
                let mut matrix: std::vec::Vec<__m128d> = [].to_vec();
                unsafe {
                    // TODO: Solve the nodal stiffness kab = \int BaDBb J
                    // Just threw in a SSE command as proof of concept.
                    matrix = vec![nodes[i], _mm_mul_pd(nodes[i], _mm_set1_pd(2.0))];
                }
                LocalStiffness {
                    matrix,
                    index: index as usize,
                }
            })
            .map(|ke| unsafe {
                // TODO: A proper global assembly
                // Just do some mm commands to test
                stiffness[ke.index] = _mm_cvtsd_f64(ke.matrix[0]);
                let temp = _mm_shuffle_pd(ke.matrix[0], ke.matrix[0], 1);
                stiffness[ke.index * 2 + 1] = _mm_cvtsd_f64(temp);
            })
            .collect::<Vec<_>>();

    // Construct force vector and map to the provided solution vector.
    // then move into memory
    // TODO: Force stream calculation and assembly like stiffness.
    // let force_stream = stream::iter_ok::<_, ()>(nodes).map().fold().map();
    let mut forces: std::vec::Vec<f64> = vec![1.0; equations];

    /* TODO: Stream stiffness calculation concurrently with force calculation.
    async {
        join!(stiffness_stream.await);
    };
    */

    // TODO: Condense this into a single function with displacement extraction.
    // Extract the forces vector
    {
        let mut i = 0;
        let mut j = 0;
        for node in 0..node_count {
            for dim in 0..2 {
                let index = (node * 2 + dim) as isize;
                if displacement_nodes.contains(&vec![node, dim]) {
                    unsafe {
                        // TODO: Set this to the calculated force.
                        *reaction_result.offset(index) = displacement_given[i];
                    }
                    i += 1;
                } else {
                    unsafe {
                        *reaction_result.offset(index) = -forces[j];
                    }
                    j += 1;
                }
            }
        }
    }

    // Do cholesky decomposition mkl::potrs
    // TODO: Debug when info indicates a failure.
    unsafe {
        let mut info = lapacke::dpotrf(
            lapacke::Layout::RowMajor,
            b'U',
            equations as i32,
            &mut stiffness.as_mut_slice(),
            equations as i32,
        );
        info = lapacke::dpotrs(
            lapacke::Layout::RowMajor,
            b'U',
            equations as i32,
            1,
            &stiffness.as_slice(),
            equations as i32,
            &mut forces.as_mut_slice(),
            1,
        );
    }

    // Calculations are in place, so forces are now displacements.
    {
        let mut i = 0;
        let mut j = 0;
        for node in 0..node_count {
            for dim in 0..2 {
                let index = (node * 2 + dim) as isize;
                if displacement_nodes.contains(&vec![node, dim]) {
                    unsafe {
                        *displacement_result.offset(index) = displacement_given[i];
                    }
                    i += 1;
                } else {
                    unsafe {
                        *displacement_result.offset(index) = -forces[j];
                    }
                    j += 1;
                }
            }
        }
    }
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
