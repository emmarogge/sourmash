pub mod minhash;
pub mod modhash;
pub mod nodegraph;

pub mod ukhs;

use serde_derive::{Deserialize, Serialize};

use crate::sketch::minhash::KmerMinHash;
use crate::sketch::modhash::KmerModHash;
use crate::sketch::ukhs::FlatUKHS;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Sketch {
    MinHash(KmerMinHash),
    ModHash(KmerModHash),
    UKHS(FlatUKHS), // FIXME
}
