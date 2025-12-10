//Program should read a genome file, substring it into reads and finally write it all to a "tirp" (tab-indexed-read-paired) file.

use std::collections::HashMap;
use std::fs;
use std::env;
use std::fs::File;
use std::fs::ReadDir;
use std::path::PathBuf;
use rand::{rng, Rng};
use rand_distr::{Normal, Distribution, Poisson};
use std::io::Write;


fn main() -> std::io::Result<()> {
    // Preamble for parse fasta
    println!("Starting program!");
    let args: Vec<String> = env::args().collect();
    let dir_genomes = std::fs::read_dir("ncbi_dataset/data").unwrap();

    // Preamble for mass_shear_seq_amp
    let min_len: usize = 250;
    let max_len: usize = 550;
    let fragm_per_bp: f64 = 0.01;
    let frag_len_distr: Normal<f32> = Normal::new(400.0, 50.0).unwrap();
    let poi = Poisson::new(2.0).unwrap();

    // Preamble for write-to-tirp
    let out_path = &args[1];
    let output_file  = std::fs::OpenOptions::new().write(true).truncate(true).create(true).open(out_path).expect("This path could NOT be used as output path. Maybe dir does not exist?");
    let phred_score_prea = (0..150).map(|_| "F").collect::<String>();
    let base_comp_prea = HashMap::from([('A', 'T'), ('T','A'),('G', 'C'), ('C','G'), ('N', 'N'),]);

    // Preamble for metafile_line_writer
    let out_metafile_path = &args[2];
    let out_metafile = std::fs::OpenOptions::new().write(true).truncate(true).create(true).open(out_metafile_path).expect("This path could NOT be used as Metafile output path. Maybe dir does not exist?");

    for (cellid_hashnum, entry) in dir_genomes.enumerate() {
      let cellid_hashnum = cellid_hashnum + 1;
      let mut assembly_dir = fs::read_dir(entry.expect("Could not read dir-entry.").path()).unwrap();
      
      let fasta_string = fs::read_to_string(find_file_in_asmbly(&mut assembly_dir, ".fna")).expect("Could not read into string");

/*       let seq_report = find_file_in_asmbly(&mut assembly_dir, "sequence_report.jsonl");
      let asmbly_report = find_file_in_asmbly(&mut assembly_dir, "assembly_data_report.jsonl"); */
    
      // Function calls
      let fna_to_hashmap  = parse_fasta(fasta_string);

      mass_seq_shear_amp(fna_to_hashmap, min_len, max_len, fragm_per_bp, frag_len_distr, &output_file, &base_comp_prea, &phred_score_prea, cellid_hashnum, &out_metafile, poi);
    }
    Ok(())
    }




fn parse_fasta(file: String) -> HashMap<String,String> {
  // Taken from Reddit (https://www.reddit.com/r/rust/comments/r5je0y/help_parsing_a_fasta_file/) lol
    let file_hashmap = file.split(">")         
      .skip(1) // Ignores first empty split - so the empty "" that happens before the first ">"
      .map(|s| {
        let mut lines = s.split("\n");
        let contig_name = lines.next().unwrap().to_string();
        let sequence_str = lines.collect::<String>();
        print!("{}", contig_name);
        (contig_name, sequence_str)
      })
      .collect::<HashMap<String,String>>();
    file_hashmap
    // Unstable sorting because it is faster 
/*     file_vec.sort_by(|a, b| a.0.to_lowercase().cmp(&b.0.to_lowercase()));
    file_vec */
}

fn mass_seq_shear_amp(hashmap_of_seqs: HashMap<String,String>, 
                    min_len: usize,
                    max_len: usize, 
                    frag_per_bp: f64, 
                    frag_len_distr: Normal<f32>,
                    output_file: &File, 
                    base_comp: &HashMap<char,char>,
                    phred_score: &String,
                    cellid_hashnum: usize,
                    out_metafile: &File,
                    poi: Poisson<f64>,) -> () {
  // Meant to take a (heading, seq) tuple and give the heading coupled with digested fragments (random substrings) stored in a Vectors with string elements.
  // Each sample should have equal amounts of sequencing depth. That SHOULD correspond to the number of times looped over a specific genome(?)
  // -- This means that different compositions of microbiome - i.e 1 E. coli vs 2 S. Typhimurim should be specified in the fasta file used for input.
  // The chromosome should also be circular.
  // Each fragment is at max 550 bp and lowest is 150 bp. This is nice as the tirp file will have no overlap "overflow" between the R1 and R2 columns since the smallest case will be R2 = reversed(R1). 
  // NOTE: Here, we are starting from after adapter trimming giving us our reads that are ONLY from the sample dna sequence.

  for (idx_in_fasta, (contig_name, sequence_str)) in hashmap_of_seqs.iter().enumerate() {
    let seq_len = sequence_str.len();

    // Decide the Copy Number via poission and some valid meta knowledge
    let copy_number = if idx_in_fasta == 0 {1 as usize} else {poi.sample(&mut rand::rng()) as usize};

    let num_of_reads = copy_number*(frag_per_bp*(sequence_str.len() as f64)).floor() as usize;
    for _ in 1..=num_of_reads {  //125000 would roughly be equal to 50 Megabits of DNA string and 1 000 000 will roughly equal 50 MB of DNA string (in UTF-8 encoding).
      let start_seq_idx: usize = rand::rng().random_range(0..=seq_len);
      // Loop to ensure fragment is within length limits
      let fragment_len: usize = loop {
        // Mean 400 by PI choice, stdev 50 which give 99,7 of reads being correct (maybe bad?)
        let fragment_len_loop: f32 = frag_len_distr.sample(&mut rng());
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
        format_and_write_to_tirp_line((&contig_name, &fragment), &output_file, &base_comp, &phred_score, &cellid_hashnum);
      } else {
        let fragment = sequence_str[start_seq_idx..(&start_seq_idx+&fragment_len)].to_string();
        format_and_write_to_tirp_line((&contig_name, &fragment), &output_file, &base_comp, &phred_score, &cellid_hashnum);
      }
      // System feedback for every 5000th fragment
/*       if fragments.len() % 5000 == 0 {
        println!("{} fragments generated.", fragments.len())
      }     */
  }
  // Write info out to the metafile
  metafile_line_writer(&out_metafile, &contig_name, &cellid_hashnum, &copy_number, &num_of_reads,);
}
                   }

fn format_and_write_to_tirp_line(tup_contigname_fragment: (&String, &String),
                                mut output_file: &File, 
                                base_comp: &HashMap<char,char>, 
                                phred_score: &String,
                                cellid_hashnum: &usize,) -> () {
  // Purpose is to use the simulated fragments, reformat them to reads and write to a .tirp file.
  // --Thinking that maybe using the identifier as the cell-id?
  // Init an output file
  // Add to main()
  // Make cell id to be used for column of index 0 in tirp
  let mut cell_id = tup_contigname_fragment.0[..11].to_string();
  cell_id.push_str("#");
  cell_id.push_str(&(format!("{:06}", cellid_hashnum)));

  let fragment = tup_contigname_fragment.1;
  // Colidx 1 to 3 is nonsense added in written string UPDATE: I see it is 1 to 2 now.
  // Make r1 and r2, cols of idx 4 and 5
  let r1 = &fragment[..150];
  // Reverses fragment and complements the 150 bases via a HashMap
  let r2 = &fragment.chars().rev().collect::<String>()[..150].chars().map(|b|{base_comp.get(&b).copied().unwrap_or('N')}).collect::<String>();
  // q1 and q2 (colidx 6 to 7) will point to phred_score and last col is a blankspace
  let out_str_line = format!("{}\t1\t1\t{}\t{}\t{}\t{}\t \n", cell_id, r1, r2, &phred_score, &phred_score);
  let _ = output_file.write_all(out_str_line.as_bytes()).expect("Problem with writing tirp data content");
    
}

fn metafile_line_writer(mut out_metafile: &File,
                        contig_name: &String,
                        cellid_hashnum: &usize,
                        copy_number: &usize,
                        num_of_reads: &usize,
//                       seq_report: &File,
                    ) -> () {
/*  cellID	strain	copyNumber	contigName	readCount
    0001	salmonella	1	main	    12
    0001	salmonella	2	plasmid1	1234
    0001	salmonella	3	plasmid2	234 */
    let strain_name = contig_name.split(" ").skip(1).next().unwrap();
    let out_str_line = format!("{}\t{}\t{}\t{}\t{}\n", format!("#{:06}", cellid_hashnum), strain_name, copy_number, contig_name, num_of_reads,);
    let _ = out_metafile.write_all(out_str_line.as_bytes()).expect("Problem with writing metadata content");
    //println!("Metaline for {} written!", contig_name)
}

fn find_file_in_asmbly(assembly_dir: &mut ReadDir, f_name_pat: &str) -> PathBuf {
    let res = assembly_dir.find(|f| f.as_ref().unwrap().file_name().into_string().unwrap().ends_with(f_name_pat)).unwrap().expect("File could not be opened").path();
    res
}