#[inline]
pub(crate) fn fast_binary_dot_product(a: u32, b: u32) -> u32 {
    (a & b).count_ones()
}
