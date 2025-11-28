//Program should read a genome file, substring it into reads and finally write it all to a "tirp" (tab-indexed-read-paired) file.

use std::fs;
use std::env;
use rand::{rng, Rng};
use rand_distr::{Normal, Distribution};
use std::io::Write;


fn main() {
    println!("Starting program!");
    let args: Vec<String> = env::args().collect();
    let file_path: &String = &args[1];
    let out_path = &args[2];
    let fasta_string: String = fs::read_to_string(file_path)
        .expect("File could not be read - maybe wrong path?");
    let vec_from_fna_parsed  = parse_fasta(fasta_string);

    println!("Vector length: {}.\nThe identifier is {}.\nThe seq is {}.", vec_from_fna_parsed.len(), vec_from_fna_parsed[0].0, vec_from_fna_parsed[0].1);

    let simulated_fragments = mass_seq_shear_amp(vec_from_fna_parsed);

    println!("First cell is from {}.\nThe first simulated fragment is {}.\nThe first (trimmed) fragment is {} bp long.", simulated_fragments[0].0, simulated_fragments[0].1[0], simulated_fragments[0].1[0].len());

    format_and_write_to_tirp(simulated_fragments, out_path);
}

// Taken from Reddit (https://www.reddit.com/r/rust/comments/r5je0y/help_parsing_a_fasta_file/) lol

fn parse_fasta(file: String) -> Vec<(String,String)> {
    let mut file_vec = file.split(">")         
      .skip(1) // Ignores first empty split - so the empty "" that happens before the first ">"
      .map(|s| {
        let mut lines = s.split("\n");
        let header = lines.next().unwrap().to_string();
        let dna = lines.collect::<String>();
        (header, dna)
      })
      .collect::<Vec<(String,String)>>();
    // Unstable sorting because it is faster 
    file_vec.sort_by(|a, b| a.0.to_lowercase().cmp(&b.0.to_lowercase()));
    file_vec
}

fn mass_seq_shear_amp(vec_of_seqs: Vec<(String,String)>) -> Vec<(String, Vec<String>)> {
  // Meant to take a (heading, seq) tuple and give the heading coupled with digested fragments (random substrings) stored in a Vectors with string elements.
  // Each sample should have equal amounts of sequencing depth. That SHOULD correspond to the number of times looped over a specific genome(?)
  // -- This means that different compositions of microbiome - i.e 1 E. coli vs 2 S. Typhimurim should be specified in the fasta file used for input.
  // The chromosome should also be circular.
  // Each fragment is at max 550 bp and lowest is 150 bp. This is nice as the tirp file will have no overlap "overflow" between the R1 and R2 columns since the smallest case will be R2 = reversed(R1). 
  // NOTE: Here, we are starting from after adapter trimming giving us our reads that are ONLY from the sample dna sequence.

  let mut out_vec= Vec::new();
  let min_len: usize = 250;
  let max_len: usize = 550;
  let fragm_per_bp: f64 = 0.01;

  for seq_tuple in vec_of_seqs {
    let identifier = seq_tuple.0;
    let sequence_str = seq_tuple.1;
    let seq_len = sequence_str.len();
    let mut fragments = Vec::new();
    let normal = Normal::new(400.0, 50.0).unwrap();

    while (fragments.len() as f64) < fragm_per_bp*(sequence_str.len() as f64) {  //125000 would roughly be equal to 50 Megabits of DNA string and 1 000 000 will roughly equal 50 MB of DNA string
      let start_seq_idx: usize = rand::rng().random_range(0..=seq_len);
      // Loop to ensure fragment is between 250 bp and 550 bp
      let fragment_len: usize = loop {
        // Mean 400 by PI choice, stdev 50 which give 99,7 of reads being correct (maybe bad?)
        let fragment_len_loop: f32 = normal.sample(&mut rng());
        let fragment_len_loop: usize = fragment_len_loop.floor() as usize;
        if min_len < fragment_len_loop && fragment_len_loop < max_len {
          break fragment_len_loop;
        }
      };
      // Ensure circular chromosome property for sequence
      if &start_seq_idx + &fragment_len > seq_len {
        let mut fragment = String::new();
        fragment.push_str(&sequence_str[start_seq_idx..]);
        fragment.push_str(&sequence_str[..(start_seq_idx+fragment_len-seq_len)]);
        fragments.push((*fragment).to_string());
      } else {
        let fragment = &sequence_str[start_seq_idx..(&start_seq_idx+&fragment_len)];
        fragments.push((*fragment).to_string());
      }
      // System feedback for every 5000th fragment
      if fragments.len() % 5000 == 0 {
        println!("{} fragments generated.", fragments.len())
      }    
      
    }
    out_vec.push((identifier, fragments));
  }
  out_vec
}

fn format_and_write_to_tirp(vec_ident_fragments: Vec<(String, Vec<String>)>, output_path: &String) -> () {
  // Purpose is to use the simulated fragments, reformat them to reads and write to a .tirp file.
  // --Thinking that maybe using the identifier as the cell-id?
  // Init an output file
  // Add to main()
  let mut out_tirp = std::fs::OpenOptions::new().write(true).truncate(true).create(true).open(output_path).expect("This path could NOT be used as output path. Maybe dir does not exist?");

  let phred_score = (0..150).map(|_| "F").collect::<String>();
  for cell in vec_ident_fragments {
    for (idx, fragment) in cell.1.iter().enumerate() {
      // Make cell id to be used for column of index 0 in tirp
      let mut cell_id = cell.0[..11].to_string();
      cell_id.push_str("#");
      cell_id.push_str(idx.to_string().trim());
      // Colidx 1 to 3 is nonsense added in written string UPDATE: I see it is 1 to 2 now.
      // Make r1 and r2, cols of idx 4 and 5
    let r1 = &fragment[..150];
    // Unsure on this one actually, how does it look from the Illumina machine???
    let r2 = &fragment.chars().rev().collect::<String>()[..150];
    // q1 and q2 (colidx 6 to 7) will point to phred_score and last col is a blankspace
    let out_str_line = format!("{}\t1\t1\t{}\t{}\t{}\t{}\t \n", cell_id, r1, r2, &phred_score, &phred_score);
    let _ = out_tirp.write_all(out_str_line.as_bytes()).expect("Problem with writing tirp data content");
    println!("{}", out_str_line)
    }
  }
}

/* Oklart om funktioner är bäst då många scopes kommer uppstå och därmed kan saker försvinna ur heap???? Counter är ju naturligtvis att se till att allt som ska vara i stacken är saker du faktiskt behöver ha med dig.

}



}
*/