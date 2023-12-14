pub mod align;
pub mod filter;
pub mod parse;
pub mod process;
pub mod reference_library;
pub mod score;
pub mod utils;

use cap::Cap;
use std::alloc;

#[global_allocator]
pub static ALLOCATOR: Cap<alloc::System> = Cap::new(alloc::System, usize::max_value());