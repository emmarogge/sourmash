use std::cmp::Ordering;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::f64::consts::PI;
use std::iter::{Iterator, Peekable};
use std::str;

use failure::Error;
use once_cell::sync::Lazy;
use serde::de::{Deserialize, Deserializer};
use serde::ser::{Serialize, SerializeStruct, Serializer};
use serde_derive::Deserialize;
use typed_builder::TypedBuilder;

use crate::_hash_murmur;
use crate::errors::SourmashError;
use crate::signature::SigsTrait;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(u32)]
pub enum HashFunctions {
    murmur64_DNA = 1,
    murmur64_protein = 2,
    murmur64_dayhoff = 3,
    murmur64_hp = 4,
}

impl std::fmt::Display for HashFunctions {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                HashFunctions::murmur64_DNA => "dna",
                HashFunctions::murmur64_protein => "protein",
                HashFunctions::murmur64_dayhoff => "dayhoff",
                HashFunctions::murmur64_hp => "hp",
            }
        )
    }
}

impl TryFrom<&str> for HashFunctions {
    type Error = Error;

    fn try_from(moltype: &str) -> Result<Self, Self::Error> {
        match moltype.to_lowercase().as_ref() {
            "dna" => Ok(HashFunctions::murmur64_DNA),
            "dayhoff" => Ok(HashFunctions::murmur64_dayhoff),
            "hp" => Ok(HashFunctions::murmur64_hp),
            "protein" => Ok(HashFunctions::murmur64_protein),
            _ => unimplemented!(),
        }
    }
}

pub fn max_hash_for_scaled(scaled: u64) -> Option<u64> {
    match scaled {
        0 => None,
        1 => Some(u64::max_value()),
        _ => Some((u64::max_value() as f64 / scaled as f64) as u64),
    }
}

pub fn scaled_for_max_hash(max_hash: u64) -> u64 {
    match max_hash {
        0 => 0,
        _ => u64::max_value() / max_hash,
    }
}

#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(Debug, Clone, PartialEq, TypedBuilder)]
pub struct KmerMinHash {
    num: u32,
    ksize: u32,

    #[builder(default_code = "HashFunctions::murmur64_DNA")]
    pub(crate) hash_function: HashFunctions,

    #[builder(default_code = "42u64")]
    seed: u64,

    #[builder(default_code = "u64::max_value()")]
    max_hash: u64,

    #[builder(default)]
    pub(crate) mins: Vec<u64>,

    #[builder(default)]
    pub(crate) abunds: Option<Vec<u64>>,
}

impl Default for KmerMinHash {
    fn default() -> KmerMinHash {
        KmerMinHash {
            num: 1000,
            ksize: 21,
            hash_function: HashFunctions::murmur64_DNA,
            seed: 42,
            max_hash: 0,
            mins: Vec::with_capacity(1000),
            abunds: None,
        }
    }
}

impl Serialize for KmerMinHash {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = match &self.abunds {
            Some(_) => 8,
            _ => 7,
        };

        let mut partial = serializer.serialize_struct("KmerMinHash", n_fields)?;
        partial.serialize_field("num", &self.num)?;
        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("seed", &self.seed)?;
        partial.serialize_field("max_hash", &self.max_hash)?;
        partial.serialize_field("mins", &self.mins)?;
        partial.serialize_field("md5sum", &self.md5sum())?;

        if let Some(abunds) = &self.abunds {
            partial.serialize_field("abundances", abunds)?;
        }

        partial.serialize_field("molecule", &self.hash_function.to_string())?;

        partial.end()
    }
}

impl<'de> Deserialize<'de> for KmerMinHash {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct TempSig {
            num: u32,
            ksize: u32,
            seed: u64,
            max_hash: u64,
            //md5sum: String,
            mins: Vec<u64>,
            abundances: Option<Vec<u64>>,
            molecule: String,
        }

        let tmpsig = TempSig::deserialize(deserializer)?;

        let num = if tmpsig.max_hash != 0 { 0 } else { tmpsig.num };
        let hash_function = match tmpsig.molecule.to_lowercase().as_ref() {
            "protein" => HashFunctions::murmur64_protein,
            "dayhoff" => HashFunctions::murmur64_dayhoff,
            "hp" => HashFunctions::murmur64_hp,
            "dna" => HashFunctions::murmur64_DNA,
            _ => unimplemented!(), // TODO: throw error here
        };

        // This shouldn't be necessary, but at some point we
        // created signatures with unordered mins =(
        let (mins, abunds) = if let Some(abunds) = tmpsig.abundances {
            let mut values: Vec<(_, _)> = tmpsig.mins.iter().zip(abunds.iter()).collect();
            values.sort();
            let mins = values.iter().map(|(v, _)| **v).collect();
            let abunds = values.iter().map(|(_, v)| **v).collect();
            (mins, Some(abunds))
        } else {
            let mut values: Vec<_> = tmpsig.mins.into_iter().collect();
            values.sort();
            (values, None)
        };

        Ok(KmerMinHash {
            num,
            ksize: tmpsig.ksize,
            seed: tmpsig.seed,
            max_hash: tmpsig.max_hash,
            mins,
            abunds,
            hash_function,
        })
    }
}

impl KmerMinHash {
    pub fn new(
        num: u32,
        ksize: u32,
        hash_function: HashFunctions,
        seed: u64,
        max_hash: u64,
        track_abundance: bool,
    ) -> KmerMinHash {
        let mins: Vec<u64>;
        let abunds: Option<Vec<u64>>;

        if num > 0 {
            mins = Vec::with_capacity(num as usize);
        } else {
            mins = Vec::with_capacity(1000);
        }

        if track_abundance {
            abunds = Some(Vec::with_capacity(mins.capacity()));
        } else {
            abunds = None
        }

        KmerMinHash {
            num,
            ksize,
            hash_function,
            seed,
            max_hash,
            mins,
            abunds,
        }
    }

    pub fn num(&self) -> u32 {
        self.num
    }

    pub fn is_protein(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_protein
    }

    fn is_dna(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_DNA
    }

    pub fn seed(&self) -> u64 {
        self.seed
    }

    pub fn max_hash(&self) -> u64 {
        self.max_hash
    }

    pub fn md5sum(&self) -> String {
        let mut md5_ctx = md5::Context::new();
        md5_ctx.consume(self.ksize().to_string());
        self.mins
            .iter()
            .for_each(|x| md5_ctx.consume(x.to_string()));
        format!("{:x}", md5_ctx.compute())
    }

    pub fn add_hash(&mut self, hash: u64) {
        self.add_hash_with_abundance(hash, 1);
    }

    pub fn add_hash_with_abundance(&mut self, hash: u64, abundance: u64) {
        let current_max = match self.mins.last() {
            Some(&x) => x,
            None => u64::max_value(),
        };

        if hash > self.max_hash && self.max_hash != 0 {
            // This is a scaled minhash, and we don't need to add the new hash
            return;
        }

        if self.num == 0 && self.max_hash == 0 {
            // why did you create this minhash? it will always be empty...
            return;
        }

        if abundance == 0 {
            // well, don't add it.
            return;
        }

        // From this point on, hash is within scaled (or no scaled specified).

        // empty mins? add it.
        if self.mins.is_empty() {
            self.mins.push(hash);
            if let Some(ref mut abunds) = self.abunds {
                abunds.push(abundance);
            }
            return;
        }

        if hash <= self.max_hash || hash <= current_max || (self.mins.len() as u32) < self.num {
            // "good" hash - within range, smaller than current entry, or
            // still have space available
            let pos = match self.mins.binary_search(&hash) {
                Ok(p) => p,
                Err(p) => p,
            };

            if pos == self.mins.len() {
                // at end - must still be growing, we know the list won't
                // get too long
                self.mins.push(hash);
                if let Some(ref mut abunds) = self.abunds {
                    abunds.push(abundance);
                }
            } else if self.mins[pos] != hash {
                // didn't find hash in mins, so inserting somewhere
                // in the middle; shrink list if needed.
                self.mins.insert(pos, hash);
                if let Some(ref mut abunds) = self.abunds {
                    abunds.insert(pos, abundance);
                }

                // is it too big now?
                if self.num != 0 && self.mins.len() > (self.num as usize) {
                    self.mins.pop();
                    if let Some(ref mut abunds) = self.abunds {
                        abunds.pop();
                    }
                }
            } else if let Some(ref mut abunds) = self.abunds {
                // pos == hash: hash value already in mins, inc count by abundance
                abunds[pos] += abundance;
            }
        }
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let hash = _hash_murmur(word, self.seed);
        self.add_hash(hash);
    }

    pub fn remove_hash(&mut self, hash: u64) {
        if let Ok(pos) = self.mins.binary_search(&hash) {
            if self.mins[pos] == hash {
                self.mins.remove(pos);
                if let Some(ref mut abunds) = self.abunds {
                    abunds.remove(pos);
                }
            }
        };
    }

    pub fn remove_many(&mut self, hashes: &[u64]) -> Result<(), Error> {
        for min in hashes {
            self.remove_hash(*min);
        }
        Ok(())
    }

    pub fn merge(&mut self, other: &KmerMinHash) -> Result<(), Error> {
        self.check_compatible(other)?;
        let max_size = self.mins.len() + other.mins.len();
        let mut merged: Vec<u64> = Vec::with_capacity(max_size);
        let mut merged_abunds: Vec<u64> = Vec::with_capacity(max_size);

        {
            let mut self_iter = self.mins.iter();
            let mut other_iter = other.mins.iter();

            let mut self_abunds_iter: Option<std::slice::Iter<'_, u64>>;
            if let Some(ref mut abunds) = self.abunds {
                self_abunds_iter = Some(abunds.iter());
            } else {
                self_abunds_iter = None;
            }

            let mut other_abunds_iter: Option<std::slice::Iter<'_, u64>>;
            if let Some(ref abunds) = other.abunds {
                other_abunds_iter = Some(abunds.iter());
            } else {
                other_abunds_iter = None;
            }

            let mut self_value = self_iter.next();
            let mut other_value = other_iter.next();
            while self_value.is_some() {
                let value = self_value.unwrap();
                match other_value {
                    None => {
                        merged.push(*value);
                        merged.extend(self_iter);
                        if let Some(sai) = self_abunds_iter {
                            merged_abunds.extend(sai);
                        }
                        break;
                    }
                    Some(x) if x < value => {
                        merged.push(*x);
                        other_value = other_iter.next();

                        if let Some(ref mut oai) = other_abunds_iter {
                            if let Some(v) = oai.next() {
                                merged_abunds.push(*v)
                            }
                        }
                    }
                    Some(x) if x == value => {
                        merged.push(*x);
                        other_value = other_iter.next();
                        self_value = self_iter.next();

                        if let Some(ref mut oai) = other_abunds_iter {
                            if let Some(v) = oai.next() {
                                if let Some(ref mut sai) = self_abunds_iter {
                                    if let Some(s) = sai.next() {
                                        merged_abunds.push(*v + *s)
                                    }
                                }
                            }
                        }
                    }
                    Some(x) if x > value => {
                        merged.push(*value);
                        self_value = self_iter.next();

                        if let Some(ref mut sai) = self_abunds_iter {
                            if let Some(v) = sai.next() {
                                merged_abunds.push(*v)
                            }
                        }
                    }
                    Some(_) => {}
                }
            }
            if let Some(value) = other_value {
                merged.push(*value);
            }
            merged.extend(other_iter);
            if let Some(oai) = other_abunds_iter {
                merged_abunds.extend(oai);
            }
        }

        if merged.len() < (self.num as usize) || (self.num as usize) == 0 {
            self.mins = merged;
            self.abunds = if merged_abunds.is_empty() {
                None
            } else {
                Some(merged_abunds)
            };
        } else {
            self.mins = merged.into_iter().take(self.num as usize).collect();
            self.abunds = if merged_abunds.is_empty() {
                None
            } else {
                Some(merged_abunds.into_iter().take(self.num as usize).collect())
            }
        }
        Ok(())
    }

    pub fn add_from(&mut self, other: &KmerMinHash) -> Result<(), Error> {
        for min in &other.mins {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many(&mut self, hashes: &[u64]) -> Result<(), Error> {
        for min in hashes {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many_with_abund(&mut self, hashes: &[(u64, u64)]) -> Result<(), Error> {
        for item in hashes {
            self.add_hash_with_abundance(item.0, item.1);
        }
        Ok(())
    }

    pub fn count_common(&self, other: &KmerMinHash, downsample: bool) -> Result<u64, Error> {
        if downsample && self.max_hash != other.max_hash {
            let (first, second) = if self.max_hash < other.max_hash {
                (self, other)
            } else {
                (other, self)
            };
            let downsampled_mh = second.downsample_max_hash(first.max_hash)?;
            first.count_common(&downsampled_mh, false)
        } else {
            self.check_compatible(other)?;
            let iter = Intersection::new(self.mins.iter(), other.mins.iter());

            Ok(iter.count() as u64)
        }
    }

    pub fn intersection(&self, other: &KmerMinHash) -> Result<(Vec<u64>, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.hash_function,
            self.seed,
            self.max_hash,
            self.abunds.is_some(),
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection::new(self.mins.iter(), other.mins.iter());

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let it2 = Intersection::new(i1.iter(), combined_mh.mins.iter());

        let common: Vec<u64> = it2.cloned().collect();
        Ok((common, combined_mh.mins.len() as u64))
    }

    pub fn intersection_size(&self, other: &KmerMinHash) -> Result<(u64, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.hash_function,
            self.seed,
            self.max_hash,
            self.abunds.is_some(),
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection::new(self.mins.iter(), other.mins.iter());

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let it2 = Intersection::new(i1.iter(), combined_mh.mins.iter());

        Ok((it2.count() as u64, combined_mh.mins.len() as u64))
    }

    // calculate Jaccard similarity, ignoring abundance.
    pub fn jaccard(&self, other: &KmerMinHash) -> Result<f64, Error> {
        self.check_compatible(other)?;
        if let Ok((common, size)) = self.intersection_size(other) {
            Ok(common as f64 / u64::max(1, size) as f64)
        } else {
            Ok(0.0)
        }
    }

    // compare two minhashes, with abundance;
    // calculate their angular similarity.
    pub fn angular_similarity(&self, other: &KmerMinHash) -> Result<f64, Error> {
        self.check_compatible(other)?;

        if self.abunds.is_none() || other.abunds.is_none() {
            // TODO: throw error, we need abundance for this
            unimplemented!() // @CTB fixme
        }

        let abunds = self.abunds.as_ref().unwrap();
        let other_abunds = other.abunds.as_ref().unwrap();

        let mut prod = 0;
        let mut other_iter = other.mins.iter().enumerate();
        let mut next_hash = other_iter.next();
        let a_sq: u64 = abunds.iter().map(|a| (a * a)).sum();
        let b_sq: u64 = other_abunds.iter().map(|a| (a * a)).sum();

        for (i, hash) in self.mins.iter().enumerate() {
            while let Some((j, k)) = next_hash {
                match k.cmp(hash) {
                    Ordering::Less => next_hash = other_iter.next(),
                    Ordering::Equal => {
                        // Calling `get_unchecked` here is safe since
                        // both `i` and `j` are valid indices
                        // (`i` and `j` came from valid iterator calls)
                        unsafe {
                            prod += abunds.get_unchecked(i) * other_abunds.get_unchecked(j);
                        }
                        break;
                    }
                    Ordering::Greater => break,
                }
            }
        }

        let norm_a = (a_sq as f64).sqrt();
        let norm_b = (b_sq as f64).sqrt();

        if norm_a == 0. || norm_b == 0. {
            return Ok(0.0);
        }
        let prod = f64::min(prod as f64 / (norm_a * norm_b), 1.);
        let distance = 2. * prod.acos() / PI;
        Ok(1. - distance)
    }

    pub fn similarity(
        &self,
        other: &KmerMinHash,
        ignore_abundance: bool,
        downsample: bool,
    ) -> Result<f64, Error> {
        if downsample && self.max_hash != other.max_hash {
            let (first, second) = if self.max_hash < other.max_hash {
                (self, other)
            } else {
                (other, self)
            };
            let downsampled_mh = second.downsample_max_hash(first.max_hash)?;
            first.similarity(&downsampled_mh, ignore_abundance, false)
        } else if ignore_abundance || self.abunds.is_none() || other.abunds.is_none() {
            self.jaccard(&other)
        } else {
            self.angular_similarity(&other)
        }
    }

    pub fn dayhoff(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_dayhoff
    }

    pub fn hp(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_hp
    }

    pub fn hash_function(&self) -> HashFunctions {
        self.hash_function
    }

    pub fn mins(&self) -> Vec<u64> {
        self.mins.clone()
    }

    // create a downsampled copy of self
    fn downsample_max_hash(&self, max_hash: u64) -> Result<KmerMinHash, Error> {
        let mut new_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.hash_function,
            self.seed,
            max_hash, // old max_hash => max_hash arg
            self.abunds.is_some(),
        );
        if self.abunds.is_some() {
            new_mh.add_many_with_abund(&self.to_vec_abunds())?;
        } else {
            new_mh.add_many(&self.mins)?;
        }
        Ok(new_mh)
    }

    fn to_vec_abunds(&self) -> Vec<(u64, u64)> {
        if let Some(abunds) = &self.abunds {
            self.mins
                .iter()
                .cloned()
                .zip(abunds.iter().cloned())
                .collect()
        } else {
            self.mins
                .iter()
                .cloned()
                .zip(std::iter::repeat(1))
                .collect()
        }
    }
}

impl SigsTrait for KmerMinHash {
    fn size(&self) -> usize {
        self.mins.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.mins.clone()
    }

    fn ksize(&self) -> usize {
        self.ksize as usize
    }

    fn check_compatible(&self, other: &KmerMinHash) -> Result<(), Error> {
        /*
        if self.num != other.num {
            return Err(SourmashError::MismatchNum {
                n1: self.num,
                n2: other.num,
            }
            .into());
        }
        */
        if self.ksize != other.ksize {
            return Err(SourmashError::MismatchKSizes.into());
        }
        if self.hash_function != other.hash_function {
            // TODO: fix this error
            return Err(SourmashError::MismatchDNAProt.into());
        }
        if self.max_hash != other.max_hash {
            return Err(SourmashError::MismatchScaled.into());
        }
        if self.seed != other.seed {
            return Err(SourmashError::MismatchSeed.into());
        }
        Ok(())
    }

    fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        let ksize = self.ksize as usize;
        let len = seq.len();

        if len < ksize {
            return Ok(());
        };

        // Here we convert the sequence to upper case and
        // pre-calculate the reverse complement for the full sequence...
        let sequence = seq.to_ascii_uppercase();
        let rc = revcomp(&sequence);

        if self.is_dna() {
            let mut last_position_check = 0;

            let mut is_valid_kmer = |i| {
                for j in std::cmp::max(i, last_position_check)..i + ksize {
                    if !VALID[sequence[j] as usize] {
                        return false;
                    }
                    last_position_check += 1;
                }
                true
            };

            for i in 0..=len - ksize {
                // ... and then while moving the k-mer window forward for the sequence
                // we move another window backwards for the RC.
                //   For a ksize = 3, and a sequence AGTCGT (len = 6):
                //                   +-+---------+---------------+-------+
                //   seq      RC     |i|i + ksize|len - ksize - i|len - i|
                //  AGTCGT   ACGACT  +-+---------+---------------+-------+
                //  +->         +->  |0|    2    |       3       |   6   |
                //   +->       +->   |1|    3    |       2       |   5   |
                //    +->     +->    |2|    4    |       1       |   4   |
                //     +->   +->     |3|    5    |       0       |   3   |
                //                   +-+---------+---------------+-------+
                // (leaving this table here because I had to draw to
                //  get the indices correctly)

                let kmer = &sequence[i..i + ksize];

                if !is_valid_kmer(i) {
                    if !force {
                        // throw error if DNA is not valid
                        return Err(SourmashError::InvalidDNA {
                            message: String::from_utf8(kmer.to_vec()).unwrap(),
                        }
                        .into());
                    }

                    continue; // skip invalid k-mer
                }

                let krc = &rc[len - ksize - i..len - i];
                self.add_word(std::cmp::min(kmer, krc));
            }
        } else {
            // protein
            let aa_ksize = self.ksize / 3;

            for i in 0..3 {
                let substr: Vec<u8> = sequence
                    .iter()
                    .cloned()
                    .skip(i)
                    .take(sequence.len() - i)
                    .collect();
                let aa = to_aa(&substr, self.dayhoff(), self.hp()).unwrap();

                aa.windows(aa_ksize as usize).for_each(|n| self.add_word(n));

                let rc_substr: Vec<u8> = rc.iter().cloned().skip(i).take(rc.len() - i).collect();
                let aa_rc = to_aa(&rc_substr, self.dayhoff(), self.hp()).unwrap();

                aa_rc
                    .windows(aa_ksize as usize)
                    .for_each(|n| self.add_word(n));
            }
        }

        Ok(())
    }

    fn add_protein(&mut self, seq: &[u8]) -> Result<(), Error> {
        let ksize = (self.ksize / 3) as usize;
        let len = seq.len();

        if len < ksize {
            return Ok(());
        }

        if let HashFunctions::murmur64_protein = self.hash_function {
            for aa_kmer in seq.windows(ksize) {
                self.add_word(&aa_kmer);
            }
            return Ok(());
        }

        let aa_seq: Vec<_> = match self.hash_function {
            HashFunctions::murmur64_dayhoff => seq.iter().cloned().map(aa_to_dayhoff).collect(),
            HashFunctions::murmur64_hp => seq.iter().cloned().map(aa_to_hp).collect(),
            invalid => {
                return Err(SourmashError::InvalidHashFunction {
                    function: format!("{}", invalid),
                }
                .into())
            }
        };

        for aa_kmer in aa_seq.windows(ksize) {
            self.add_word(aa_kmer);
        }

        Ok(())
    }
}

struct Intersection<T, I: Iterator<Item = T>> {
    iter: Peekable<I>,
    other: Peekable<I>,
}

impl<T, I: Iterator<Item = T>> Intersection<T, I> {
    pub fn new(left: I, right: I) -> Self {
        Intersection {
            iter: left.peekable(),
            other: right.peekable(),
        }
    }
}

impl<T: Ord, I: Iterator<Item = T>> Iterator for Intersection<T, I> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            let res = match (self.iter.peek(), self.other.peek()) {
                (Some(ref left_key), Some(ref right_key)) => left_key.cmp(right_key),
                _ => return None,
            };

            match res {
                Ordering::Less => {
                    self.iter.next();
                }
                Ordering::Greater => {
                    self.other.next();
                }
                Ordering::Equal => {
                    self.other.next();
                    return self.iter.next();
                }
            }
        }
    }
}

const COMPLEMENT: [u8; 256] = {
    let mut lookup = [0; 256];
    lookup[b'A' as usize] = b'T';
    lookup[b'C' as usize] = b'G';
    lookup[b'G' as usize] = b'C';
    lookup[b'T' as usize] = b'A';
    lookup[b'N' as usize] = b'N';
    lookup
};

#[inline]
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|nt| COMPLEMENT[*nt as usize])
        .collect()
}

static CODONTABLE: Lazy<HashMap<&'static str, u8>> = Lazy::new(|| {
    [
        // F
        ("TTT", b'F'),
        ("TTC", b'F'),
        // L
        ("TTA", b'L'),
        ("TTG", b'L'),
        // S
        ("TCT", b'S'),
        ("TCC", b'S'),
        ("TCA", b'S'),
        ("TCG", b'S'),
        ("TCN", b'S'),
        // Y
        ("TAT", b'Y'),
        ("TAC", b'Y'),
        // *
        ("TAA", b'*'),
        ("TAG", b'*'),
        // *
        ("TGA", b'*'),
        // C
        ("TGT", b'C'),
        ("TGC", b'C'),
        // W
        ("TGG", b'W'),
        // L
        ("CTT", b'L'),
        ("CTC", b'L'),
        ("CTA", b'L'),
        ("CTG", b'L'),
        ("CTN", b'L'),
        // P
        ("CCT", b'P'),
        ("CCC", b'P'),
        ("CCA", b'P'),
        ("CCG", b'P'),
        ("CCN", b'P'),
        // H
        ("CAT", b'H'),
        ("CAC", b'H'),
        // Q
        ("CAA", b'Q'),
        ("CAG", b'Q'),
        // R
        ("CGT", b'R'),
        ("CGC", b'R'),
        ("CGA", b'R'),
        ("CGG", b'R'),
        ("CGN", b'R'),
        // I
        ("ATT", b'I'),
        ("ATC", b'I'),
        ("ATA", b'I'),
        // M
        ("ATG", b'M'),
        // T
        ("ACT", b'T'),
        ("ACC", b'T'),
        ("ACA", b'T'),
        ("ACG", b'T'),
        ("ACN", b'T'),
        // N
        ("AAT", b'N'),
        ("AAC", b'N'),
        // K
        ("AAA", b'K'),
        ("AAG", b'K'),
        // S
        ("AGT", b'S'),
        ("AGC", b'S'),
        // R
        ("AGA", b'R'),
        ("AGG", b'R'),
        // V
        ("GTT", b'V'),
        ("GTC", b'V'),
        ("GTA", b'V'),
        ("GTG", b'V'),
        ("GTN", b'V'),
        // A
        ("GCT", b'A'),
        ("GCC", b'A'),
        ("GCA", b'A'),
        ("GCG", b'A'),
        ("GCN", b'A'),
        // D
        ("GAT", b'D'),
        ("GAC", b'D'),
        // E
        ("GAA", b'E'),
        ("GAG", b'E'),
        // G
        ("GGT", b'G'),
        ("GGC", b'G'),
        ("GGA", b'G'),
        ("GGG", b'G'),
        ("GGN", b'G'),
    ]
    .iter()
    .cloned()
    .collect()
});

// Dayhoff table from
// Peris, P., López, D., & Campos, M. (2008).
// IgTM: An algorithm to predict transmembrane domains and topology in
// proteins. BMC Bioinformatics, 9(1), 1029–11.
// http://doi.org/10.1186/1471-2105-9-367
//
// Original source:
// Dayhoff M. O., Schwartz R. M., Orcutt B. C. (1978).
// A model of evolutionary change in proteins,
// in Atlas of Protein Sequence and Structure,
// ed Dayhoff M. O., editor.
// (Washington, DC: National Biomedical Research Foundation; ), 345–352.
//
// | Amino acid    | Property              | Dayhoff |
// |---------------|-----------------------|---------|
// | C             | Sulfur polymerization | a       |
// | A, G, P, S, T | Small                 | b       |
// | D, E, N, Q    | Acid and amide        | c       |
// | H, K, R       | Basic                 | d       |
// | I, L, M, V    | Hydrophobic           | e       |
// | F, W, Y       | Aromatic              | f       |
static DAYHOFFTABLE: Lazy<HashMap<u8, u8>> = Lazy::new(|| {
    [
        // a
        (b'C', b'a'),
        // b
        (b'A', b'b'),
        (b'G', b'b'),
        (b'P', b'b'),
        (b'S', b'b'),
        (b'T', b'b'),
        // c
        (b'D', b'c'),
        (b'E', b'c'),
        (b'N', b'c'),
        (b'Q', b'c'),
        // d
        (b'H', b'd'),
        (b'K', b'd'),
        (b'R', b'd'),
        // e
        (b'I', b'e'),
        (b'L', b'e'),
        (b'M', b'e'),
        (b'V', b'e'),
        // e
        (b'F', b'f'),
        (b'W', b'f'),
        (b'Y', b'f'),
    ]
    .iter()
    .cloned()
    .collect()
});

// HP Hydrophobic/hydrophilic mapping
// From: Phillips, R., Kondev, J., Theriot, J. (2008).
// Physical Biology of the Cell. New York: Garland Science, Taylor & Francis Group. ISBN: 978-0815341635

//
// | Amino acid                            | HP
// |---------------------------------------|---------|
// | A, F, G, I, L, M, P, V, W, Y          | h       |
// | N, C, S, T, D, E, R, H, K, Q          | p       |
static HPTABLE: Lazy<HashMap<u8, u8>> = Lazy::new(|| {
    [
        // h
        (b'A', b'h'),
        (b'F', b'h'),
        (b'G', b'h'),
        (b'I', b'h'),
        (b'L', b'h'),
        (b'M', b'h'),
        (b'P', b'h'),
        (b'V', b'h'),
        (b'W', b'h'),
        (b'Y', b'h'),
        // p
        (b'N', b'p'),
        (b'C', b'p'),
        (b'S', b'p'),
        (b'T', b'p'),
        (b'D', b'p'),
        (b'E', b'p'),
        (b'R', b'p'),
        (b'H', b'p'),
        (b'K', b'p'),
        (b'Q', b'p'),
    ]
    .iter()
    .cloned()
    .collect()
});

#[inline]
pub(crate) fn translate_codon(codon: &[u8]) -> Result<u8, Error> {
    if codon.len() == 1 {
        return Ok(b'X');
    }

    if codon.len() == 2 {
        let mut v = codon.to_vec();
        v.push(b'N');
        match CODONTABLE.get(str::from_utf8(v.as_slice()).unwrap()) {
            Some(aa) => return Ok(*aa),
            None => return Ok(b'X'),
        }
    }

    if codon.len() == 3 {
        match CODONTABLE.get(str::from_utf8(codon).unwrap()) {
            Some(aa) => return Ok(*aa),
            None => return Ok(b'X'),
        }
    }

    Err(SourmashError::InvalidCodonLength {
        message: format!("{}", codon.len()),
    }
    .into())
}

#[inline]
pub(crate) fn aa_to_dayhoff(aa: u8) -> u8 {
    match DAYHOFFTABLE.get(&aa) {
        Some(letter) => *letter,
        None => b'X',
    }
}

pub(crate) fn aa_to_hp(aa: u8) -> u8 {
    match HPTABLE.get(&aa) {
        Some(letter) => *letter,
        None => b'X',
    }
}

#[inline]
fn to_aa(seq: &[u8], dayhoff: bool, hp: bool) -> Result<Vec<u8>, Error> {
    let mut converted: Vec<u8> = Vec::with_capacity(seq.len() / 3);

    for chunk in seq.chunks(3) {
        if chunk.len() < 3 {
            break;
        }

        let residue = translate_codon(chunk)?;
        if dayhoff {
            converted.push(aa_to_dayhoff(residue) as u8);
        } else if hp {
            converted.push(aa_to_hp(residue) as u8);
        } else {
            converted.push(residue);
        }
    }

    Ok(converted)
}

const VALID: [bool; 256] = {
    let mut lookup = [false; 256];
    lookup[b'A' as usize] = true;
    lookup[b'C' as usize] = true;
    lookup[b'G' as usize] = true;
    lookup[b'T' as usize] = true;
    lookup
};
