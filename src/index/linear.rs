use std::path::Path;
use std::rc::Rc;

use derive_builder::Builder;
use failure::Error;

use crate::index::storage::Storage;
use crate::index::{Comparable, Index};

#[derive(Builder)]
pub struct LinearIndex<L> {
    //#[builder(setter(skip))]
    storage: Rc<dyn Storage>,

    #[builder(setter(skip))]
    pub(crate) datasets: Vec<L>,
}

impl<L> LinearIndex<L>
where
    L: Default,
{
    pub fn builder() -> LinearIndexBuilder<L> {
        LinearIndexBuilder::default()
    }
}

impl<L> Index for LinearIndex<L>
where
    L: Clone + Comparable<L>,
{
    type Item = L;

    fn find<F>(
        &self,
        search_fn: F,
        sig: &Self::Item,
        threshold: f64,
    ) -> Result<Vec<&Self::Item>, Error>
    where
        F: Fn(&dyn Comparable<Self::Item>, &Self::Item, f64) -> bool,
    {
        Ok(self
            .datasets
            .iter()
            .flat_map(|node| {
                if search_fn(node, sig, threshold) {
                    Some(node)
                } else {
                    None
                }
            })
            .collect())
    }

    fn insert(&mut self, node: &L) -> Result<(), Error> {
        self.datasets.push(node.clone());
        Ok(())
    }

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        Ok(())
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        Ok(())
    }

    fn datasets(&self) -> Vec<Self::Item> {
        self.datasets.to_vec()
    }
}
