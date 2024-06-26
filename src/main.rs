use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{stdin, BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

fn main() {
    let file = env::args()
        .nth(1)
        .expect("Please enter a existing file path (arg 1)");
    if file == "-help" {
        println!(
            "
Mokuba v1.0.0

Description: This takes a fasta file and searches by the id to find the corresponding sequence.

Format: mokuba /path/to/fasta search -flags

Example: mokuba medtr.A17.gnm5.ann1_6.L2RX.cds.fna.gz Chr1g0147651 -md

Flags:
    -s: For fastas with a single sequence
    -sd: For fastas with a single sequence that need to be unziped from .gz
    -m: For fastas with multiple sequences
    -md: For fastas with multiple sequences that need to be unziped from .gz
"
        );
        std::process::exit(3);
    }

    let id = env::args()
        .nth(2)
        .expect("Please enter an id to search for (arg 2)");
    let option = env::args().nth(3).expect("Please enter an option (arg 3)");

    let opt_id1 = String::from("-s");
    let opt_id2 = String::from("-sd");
    let opt_id3 = String::from("-m");
    let opt_id4 = String::from("-md");

    let mut write_id = String::new();
    let mut write_seq = String::new();

    if option == opt_id1 {
        let fasta_read = read_single_fasta(file.clone());
        let info = get_info(fasta_read, id.clone());
        write_id.push_str(&info.0);
        write_seq.push_str(&info.1);
        promts(write_id.clone(), write_seq.clone());
    }
    if option == opt_id2 {
        let fasta_read = read_single_fasta_deco(file.clone());
        let info = get_info(fasta_read, id.clone());
        write_id.push_str(&info.0);
        write_seq.push_str(&info.1);
        promts(write_id.clone(), write_seq.clone());
    }
    if option == opt_id3 {
        let fasta_read = read_multiple_fasta(file.clone());
        let info = get_info(fasta_read, id.clone());
        write_id.push_str(&info.0);
        write_seq.push_str(&info.1);
        promts(write_id.clone(), write_seq.clone());
    }
    if option == opt_id4 {
        let fasta_read = read_multiple_fasta_deco(file.clone());
        let info = get_info(fasta_read, id.clone());
        write_id.push_str(&info.0);
        write_seq.push_str(&info.1);
        promts(write_id.clone(), write_seq.clone());
    }
}

fn promts(write_id: String, write_seq: String) {
    let mut write_file = String::new();
    println!("\nWould you like to write the output to a file? (Y/N)\n");
    stdin()
        .read_line(&mut write_file)
        .expect("Could not read entry");
    let write_file = write_file.trim();
    if write_file.to_uppercase() == "Y" {
        let mut file_name = String::new();
        println!("\nEnter a name for the file\n");
        stdin()
            .read_line(&mut file_name)
            .expect("Could not read entry");
        let file_name = file_name.trim();
        write_seq_file(file_name, write_id, write_seq);
    }
    println!("\nDone\n");
}

fn write_seq_file(name: &str, id: String, seq: String) {
    let mut file_name = String::from(name);
    file_name.push_str(".fasta");
    if Path::new(&file_name).exists() {
        println!(
            "{} already exists, please move it to create new meta data file",
            &file_name
        );
    } else {
        File::create(&file_name).expect("Could not create file");
        let mut file_write = File::options().append(true).open(&file_name).unwrap();
        writeln!(&mut file_write, "{}", id).expect("Could not write ID");
        let mut counter = 0;
        let seq_chars = seq.chars();
        let mut seq_string = String::new();
        for ch in seq_chars {
            if counter < 60 {
                counter += 1;
                seq_string.push(ch);
            } else {
                writeln!(&mut file_write, "{}", seq_string).expect("Could not write sequence");
                counter = 0;
                seq_string.clear();
            }
        }
    }
}

fn get_info(hash: HashMap<String, String>, id: String) -> (String, String) {
    let found_id = get_id(hash.clone(), id);
    let retriev_info = hash.get_key_value(&found_id).expect("Could not find id");
    let seq = retriev_info.1;
    let validate = retriev_info.0;
    assert!(&found_id == validate);
    println!("\nid:\n{:?}\n", found_id);
    println!("sequence:\n{:?}\n", seq);
    (found_id, seq.to_string())
}

fn get_id(hash: HashMap<String, String>, search_id: String) -> String {
    let mut found_id = String::new();
    let ids = hash.keys();
    for i in ids {
        if i.contains(&search_id) {
            found_id.push_str(&i);
        }
    }
    found_id
}

fn read_single_fasta<P>(filename: P) -> HashMap<String, String>
where
    P: AsRef<Path>,
{
    let file = File::open(filename).expect("Could not open file");
    let buf = BufReader::new(file);
    let mut fasta = HashMap::new();
    let mut curid = String::new();
    let mut curseq = String::new();
    for line in buf.lines() {
        let line = line.expect("Could not read line");
        if line.starts_with('>') {
            curid = line[..].trim().to_string();
        } else {
            curseq.push_str(line.trim());
        }
    }
    fasta.insert(curid.clone(), curseq.clone());
    fasta
}

fn read_single_fasta_deco<P>(filename: P) -> HashMap<String, String>
where
    P: AsRef<Path>,
{
    let file = File::open(filename).expect("Could not open file");
    let gz = MultiGzDecoder::new(file);
    let buf = BufReader::new(gz);
    let mut fasta = HashMap::new();
    let mut curid = String::new();
    let mut curseq = String::new();
    for line in buf.lines() {
        let line = line.expect("Could not read line");
        if line.starts_with('>') {
            curid = line[..].trim().to_string();
        } else {
            curseq.push_str(line.trim());
        }
    }
    fasta.insert(curid.clone(), curseq.clone());
    fasta
}

fn read_multiple_fasta_deco<P>(filename: P) -> HashMap<String, String>
where
    P: AsRef<Path>,
{
    let start = Instant::now();
    let file = File::open(filename).expect("Could not open file");
    let gz = MultiGzDecoder::new(file);
    let buf = BufReader::new(gz);
    let mut fasta = HashMap::new();
    let mut curid = String::new();
    let mut curseq = String::new();
    for line in buf.lines() {
        let line = line.expect("Could not read line");
        if line.starts_with('>') {
            if !curid.is_empty() {
                fasta.insert(curid.clone(), curseq.clone());
                curseq.clear();
            }
            // println!("{:?}", &line[..].trim());
            curid = line[..].trim().to_string();
        } else {
            curseq.push_str(line.trim());
        }
    }
    let duration = start.elapsed();
    println!("It took {:?} to decode and read", duration);
    fasta
}

fn read_multiple_fasta<P>(filename: P) -> HashMap<String, String>
where
    P: AsRef<Path>,
{
    let start = Instant::now();
    let file = File::open(filename).expect("Could not open file");
    let buf = BufReader::new(file);
    let mut fasta = HashMap::new();
    let mut curid = String::new();
    let mut curseq = String::new();
    for line in buf.lines() {
        let line = line.expect("Could not read line");
        if line.starts_with('>') {
            if !curid.is_empty() {
                fasta.insert(curid.clone(), curseq.clone());
                curseq.clear();
            }
            // println!("{:?}", &line[..].trim());
            curid = line[..].trim().to_string();
        } else {
            curseq.push_str(line.trim());
        }
    }
    let duration = start.elapsed();
    println!("It took {:?} to decode and read", duration);
    fasta
}
