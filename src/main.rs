use flate2::read::MultiGzDecoder;
use indicatif::ProgressBar;
use std::collections::HashMap;
use std::env::{args, set_var};
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, Write};
use std::path::Path;

fn main() {
    set_var("RUST_BACKTRACE", "1");
    let option = args().nth(1).expect("Please enter an option (arg 2)");
    if option == "--help" {
        println!(
            "
Mokuba v1.0.0

Description: This takes a fasta file and searches by the id to find the corresponding sequence.

Format: mokuba -option /path/to/fasta search 

Example: mokuba md- medtr.A17.gnm5.ann1_6.L2RX.cds.fna.gz Chr1g0147651 

Flags:
    -m: For fastas with multiple sequences
    -md: For fastas with multiple sequences that need to be unziped from .gz
    -si: Does not save file, instead only gives standard out. Good for piping. Only for .gz files at the moment.
"
        );
        std::process::exit(3);
    }

    let file = args()
        .nth(2)
        .expect("Please enter a existing file path (arg 1)");

    let mut id = String::new();
    if option == "-sio" {
        let st_in = stdin().lines();
        for line in st_in {
            id.push_str(&line.unwrap());
            id.push_str(",");
        }
    } else {
        id.push_str(
            &args()
                .nth(3)
                .expect("Please enter an id to search for (arg 3)"),
        );
    }
    let opt_id1 = String::from("-m");
    let opt_id2 = String::from("-md");
    let opt_id3 = String::from("-sio");

    let mut write_id = String::new();
    let mut write_seq = String::new();

    if option == opt_id1 {
        let fasta_read = read_multiple_fasta(file.clone());
        let info = get_info(&fasta_read, &id, &option);
        write_id.push_str(&info.0);
        write_seq.push_str(&info.1);
        promts(write_id.clone(), write_seq.clone());
    }
    if option == opt_id2 {
        let fasta_read = read_multiple_fasta_deco(file.clone());
        let info = get_info(&fasta_read, &id, &option);
        write_id.push_str(&info.0);
        write_seq.push_str(&info.1);
        promts(write_id.clone(), write_seq.clone());
    }
    if option == opt_id3 {
        let fasta_read = read_multiple_fasta_deco(file.clone());
        let id_iter = id.split(",");
        let id_vec: Vec<String> = id_iter.clone().map(|i| i.to_string()).collect();
        println!("Getting sequences from ID's");
        let bar = ProgressBar::new(id_vec.len() as u64);
        for id in id_iter {
            bar.inc(1);
            get_info(&fasta_read, &id, &option);
        }
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

fn get_info(hash: &HashMap<String, String>, id: &str, option: &str) -> (String, String) {
    let found_id = get_id(hash.clone(), id);
    let retriev_info_op = hash.get_key_value(&found_id);
    let retriev_info = match retriev_info_op {
        Some(info) => info,
        None => (&String::from("NA"), &String::from("NA")),
    };

    let seq = retriev_info.1;
    if option == "-sio" {
        let mut write_id = String::from(">");
        write_id.push_str(id);
        stdout()
            .write_all(format!("{}\n{}\n", seq, write_id).as_bytes())
            .unwrap();
        (found_id, seq.to_string())
    } else {
        println!("\nid:\n");
        stdout()
            .write_all(format!("\n{}", found_id).as_bytes())
            .unwrap();
        println!("sequence:\n");
        stdout()
            .write_all(format!("\n{}\n", seq).as_bytes())
            .unwrap();
        (found_id, seq.to_string())
    }
}

fn get_id(hash: HashMap<String, String>, search_id: &str) -> String {
    let mut found_id = String::new();
    let ids = hash.keys();
    for i in ids {
        if i.contains(&search_id) {
            found_id.push_str(&i);
        }
    }
    found_id
}

fn read_multiple_fasta_deco<P>(filename: P) -> HashMap<String, String>
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
            if !curid.is_empty() {
                fasta.insert(curid.clone(), curseq.clone());
                curseq.clear();
            }
            curid = line[..].trim().to_string();
        } else {
            curseq.push_str(line.trim());
        }
    }
    fasta
}

fn read_multiple_fasta<P>(filename: P) -> HashMap<String, String>
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
            if !curid.is_empty() {
                fasta.insert(curid.clone(), curseq.clone());
                curseq.clear();
            }
            curid = line[..].trim().to_string();
        } else {
            curseq.push_str(line.trim());
        }
    }
    fasta
}
