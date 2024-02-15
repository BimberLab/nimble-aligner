pub mod align;
pub mod filter;
pub mod parse;
pub mod process;
pub mod reference_library;
pub mod score;
pub mod utils;

use tikv_jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;
