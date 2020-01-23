pub mod minhash;
pub mod nodegraph;

use std::any::Any;

use dyn_clone::DynClone;
use failure::Error;

use crate::sketch::minhash::HashFunctions;

//https://github.com/dtolnay/typetag
//https://github.com/dtolnay/typetag/issues/10
#[typetag::serde(tag = "type")]
pub trait Sketch: DynClone {
    fn size(&self) -> usize;
    fn to_vec(&self) -> Vec<u64>;
    fn check_compatible(&self, other: &Box<dyn Sketch>) -> Result<(), Error>;
    fn add_sequence(&mut self, seq: &[u8], _force: bool) -> Result<(), Error>;
    fn add_protein(&mut self, seq: &[u8]) -> Result<(), Error>;
    fn ksize(&self) -> usize;
    fn hash_function(&self) -> HashFunctions;

    fn similarity(&self, other: &Box<dyn Sketch>) -> Result<f64, Error>;
    fn containment(&self, other: &Box<dyn Sketch>) -> Result<f64, Error>;
    fn md5sum(&self) -> String;

    fn as_any(&self) -> &dyn Any;
}

dyn_clone::clone_trait_object!(Sketch);

impl std::fmt::Debug for Box<dyn Sketch> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //write!(f, "{{ x: {}, y: {} }}", self.x, self.y)
        write!(f, "Sketch")
    }
}

impl PartialEq for Box<dyn Sketch> {
    fn eq(&self, other: &Box<dyn Sketch>) -> bool {
        unimplemented!()
    }
}
