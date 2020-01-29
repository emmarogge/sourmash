use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::HashFunctions;
use sourmash::sketch::modhash::KmerModHash;

#[test]
fn throws_error() {
    let mut mh = KmerModHash::new(4, HashFunctions::murmur64_DNA, 42, 1, false);

    match mh.add_sequence(b"ATGR", false) {
        Ok(_) => assert!(false, "R is not a valid DNA character"),
        Err(_) => assert!(true),
    }
}

#[test]
fn merge() {
    let mut a = KmerModHash::new(10, HashFunctions::murmur64_DNA, 42, 2, false);
    let mut b = KmerModHash::new(10, HashFunctions::murmur64_DNA, 42, 2, false);

    a.add_sequence(b"TGCCGCCCAGCA", false).unwrap();
    b.add_sequence(b"TGCCGCCCAGCA", false).unwrap();

    a.add_sequence(b"GTCCGCCCAGTGA", false).unwrap();
    b.add_sequence(b"GTCCGCCCAGTGG", false).unwrap();

    a.merge(&b).unwrap();
    assert_eq!(
        a.to_vec(),
        vec![
            8373222269469409550,
            11085758717695534616,
            11760449009842383350,
        ]
    );
}

#[test]
fn compare() {
    let mut a = KmerModHash::new(10, HashFunctions::murmur64_DNA, 42, 21, false);
    let mut b = KmerModHash::new(10, HashFunctions::murmur64_DNA, 42, 21, false);

    a.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false)
        .unwrap();
    b.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false)
        .unwrap();
    assert_eq!(a.compare(&b).unwrap(), 1.0);
    assert_eq!(b.compare(&a).unwrap(), 1.0);

    b.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false)
        .unwrap();
    assert_eq!(a.compare(&b).unwrap(), 1.0);
    assert_eq!(b.compare(&a).unwrap(), 1.0);

    b.add_sequence(b"GATTGGTGCACACTTAACTGGGTGCCGCGCTGGTGCTGATCCATGAAGTT", false)
        .unwrap();
    assert!(a.compare(&b).unwrap() >= 0.3);
    assert!(b.compare(&a).unwrap() >= 0.3);
}

#[test]
fn invalid_dna() {
    let mut a = KmerModHash::new(3, HashFunctions::murmur64_DNA, 42, 1, false);

    a.add_sequence(b"AAANNCCCTN", true).unwrap();
    assert_eq!(a.mins().len(), 3);

    let mut b = KmerModHash::new(3, HashFunctions::murmur64_DNA, 42, 1, false);
    b.add_sequence(b"NAAA", true).unwrap();
    assert_eq!(b.mins().len(), 1);
}

#[test]
fn similarity() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = KmerModHash::new(20, HashFunctions::murmur64_hp, 42, 1, true);
    let mut b = KmerModHash::new(20, HashFunctions::murmur64_hp, 42, 1, true);

    a.add_hash(1);
    b.add_hash(1);
    b.add_hash(2);

    assert!((a.similarity(&a, false)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, false)? - 0.5).abs() < 0.001);

    Ok(())
}

#[test]
fn similarity_2() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = KmerModHash::new(5, HashFunctions::murmur64_DNA, 42, 2, true);
    let mut b = KmerModHash::new(5, HashFunctions::murmur64_DNA, 42, 2, true);

    a.add_sequence(b"ATGGA", false)?;
    a.add_sequence(b"GGACA", false)?;

    a.add_sequence(b"ATGGA", false)?;
    b.add_sequence(b"ATGGA", false)?;

    assert!(
        (a.similarity(&b, false)? - 0.705).abs() < 0.001,
        "{}",
        a.similarity(&b, false)?
    );

    Ok(())
}

#[test]
fn similarity_3() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = KmerModHash::new(20, HashFunctions::murmur64_dayhoff, 42, 2, true);
    let mut b = KmerModHash::new(20, HashFunctions::murmur64_dayhoff, 42, 2, true);

    a.add_hash(2);
    a.add_hash(2);
    a.add_hash(4);
    a.add_hash(4);

    b.add_hash(1);
    b.add_hash(2);
    b.add_hash(4);
    b.add_hash(6);

    assert!((a.similarity(&a, false)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, false)? - 0.608).abs() < 0.001);

    assert!((a.similarity(&a, true)? - 1.0).abs() < 0.001);
    assert!((a.similarity(&b, true)? - 0.66666).abs() < 0.001);

    Ok(())
}

#[test]
fn dayhoff() {
    let mut a = KmerModHash::new(6, HashFunctions::murmur64_dayhoff, 42, 1, false);
    let mut b = KmerModHash::new(6, HashFunctions::murmur64_protein, 42, 1, false);

    a.add_sequence(b"ACTGAC", false).unwrap();
    b.add_sequence(b"ACTGAC", false).unwrap();

    assert_eq!(a.size(), 2);
    assert_eq!(b.size(), 2);
}

#[test]
fn hp() {
    let mut a = KmerModHash::new(6, HashFunctions::murmur64_hp, 42, 1, false);
    let mut b = KmerModHash::new(6, HashFunctions::murmur64_protein, 42, 1, false);

    a.add_sequence(b"ACTGAC", false).unwrap();
    b.add_sequence(b"ACTGAC", false).unwrap();

    assert_eq!(a.size(), 2);
    assert_eq!(b.size(), 2);
}
