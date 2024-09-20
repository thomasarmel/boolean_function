//! g and f are affine equivalent if forall x, g(x) = f(Dx + a) + bx + c for some D n*n matrix, a, b, c vectors

use crate::SmallBooleanFunction;

pub const BOOLEAN_FUNCTIONS_3_VAR_AFFINE_EQ_CLASSES: [SmallBooleanFunction; 3] = [
    SmallBooleanFunction::from_truth_table_unchecked(0xaa, 3),
    SmallBooleanFunction::from_truth_table_unchecked(0xab, 3),
    SmallBooleanFunction::from_truth_table_unchecked(0xac, 3)
];

pub const BOOLEAN_FUNCTIONS_4_VAR_AFFINE_EQ_CLASSES: [SmallBooleanFunction; 8] = [
    SmallBooleanFunction::from_truth_table_unchecked(0xaa55, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xab55, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xbb55, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xaba5, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xaaff, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xaba4, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xab12, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xac90, 4), // bent
];

pub const BOOLEAN_FUNCTIONS_5_VAR_AFFINE_EQ_CLASSES: [SmallBooleanFunction; 48] = [
    SmallBooleanFunction::from_truth_table_unchecked(0xaa55aa55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaa55ab55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaa55bb55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaa5dbb55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaaddbb55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaa5dbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x2a5dbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaaddbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x2a5dbf51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x6a5dbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x2addbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xa8ddbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaeddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x0a5dbf51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8addda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xa8dd9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88ddbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88ddbb11, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8c5dda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xa89d9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8eddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaefdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x025dbf51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88ddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88dd9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xceddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x0eddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x425dbf51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8cddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88dddb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x289d9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x86fdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88dddb71, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xcefdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x0efdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x288d9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8cfdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8cdddb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8ccdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x289d9b41, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x488ddb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xccfdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x688d9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x288d9b41, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x288d1b41, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xdcfdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x68ad9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x688ddb51, 5),
];