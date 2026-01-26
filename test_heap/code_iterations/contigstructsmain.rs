//Program should read a genome file, substring it into reads and finally write it all to a "tirp" (tab-indexed-read-paired) file.

use std::collections::HashMap;
use std::fs;
use std::env;
use std::fs::File;
use std::fs::ReadDir;
use std::time::Duration;
use rand::Rng;
use rand_distr::{Normal, Distribution, Poisson};
use std::io::{BufWriter,Write};
use json;
use std::fs::DirEntry;
use indicatif::ProgressBar;


struct Bactdatafromfasta <'a>{
  fasta_as_vec_hashmap: Vec<(String,HashMap<&'a str,String>)>,
  strain_name: String,
  assembly_name: String,
}

struct Contigfromfasta <'a>{
  name: String,
  seq_str: String,
  chr_name: String,
  strain_name: String,
  assembly_name: String,
}

struct CellidHashnum {
    id_counter: usize
}

impl CellidHashnum {
    fn count(&mut self) -> () {
        self.id_counter = self.id_counter + 1;
    }
}

const BASE_COMP_LUT: [u8; 256] = {
      let mut lut = [0u8; 256];
      lut[b'A' as usize] = b'T';
      lut[b'T' as usize] = b'A';
      lut[b'C' as usize] = b'G';
      lut[b'G' as usize] = b'C';
      lut[b'N' as usize] = b'N';
      lut
    };

fn main() -> std::io::Result<()> {
    // Preamble for parse fasta
    println!("Starting program!");
    let args: Vec<String> = env::args().collect();
    let dir_genomes = std::fs::read_dir("ncbi_dataset/data").unwrap();

    // Preamble for read_simulator
    let frag_len_distr: Normal<f32> = Normal::new(400.0, 50.0).unwrap();
    let poi = Poisson::new(13.0).unwrap();

    // Preamble for write-to-tirp
    let mut output_file  = BufWriter::new(std::fs::OpenOptions::new().write(true).truncate(true).create(true).open(&args[1]).expect("This path could NOT be used as output path. Maybe dir does not exist?"));
    let phred_score_prea = (0..150).map(|_| "F").collect::<String>();

    // Preamble for write-chrom-to-tirp
    let bytes_chrom_r1 = *&args[3].parse::<usize>().unwrap();


    // Preamble for metafile_line_writer
    let mut out_metafile = BufWriter::new(std::fs::OpenOptions::new().write(true).truncate(true).create(true).open(&args[2]).expect("This path could NOT be used as Metafile output path. Maybe dir does not exist?"));

    let metagenome_json  = json::parse(&fs::read_to_string("run-setup.json").expect("could not read run-setup.json!")).expect("Problem with run-setup.json");

    let mut cellid_hashnum = CellidHashnum {
      id_counter: 1,
    };

    for entry in dir_genomes {
      let bact_entry = instantiate_bact(entry.expect("Problem with reading entry!"));
      // Function (all things for an entry so it is loaded and ready) -> <HashMap(PathBuf, HashMap(String, ))
      let num_bacteria = (&metagenome_json[&bact_entry.assembly_name].to_string()).parse::<usize>().expect("Could not convert into usize.");
      for _ in 0..num_bacteria {
        // Function calls
        read_simulator(&bact_entry, frag_len_distr, &mut output_file, &phred_score_prea, cellid_hashnum.id_counter, &mut out_metafile, poi, bytes_chrom_r1);
        cellid_hashnum.count()
      }
    }
    output_file.flush().expect("Problem with flushing output tirp file!");
    out_metafile.flush().expect("Problem with flushing output metafile!");
    Ok(())
    }

// Simulation function
fn read_simulator(bact_entry: &Bactdatafromfasta,
                    frag_len_distr: Normal<f32>,
                    output_file: &mut BufWriter<File>,
                    phred_score: &String,
                    cellid_hashnum: usize,
                    out_metafile: &mut BufWriter<File>,
                    poi: Poisson<f64>,
                    bytes_chrom_r1: usize
                ) -> () {
  // Meant to take a (heading, seq) tuple and give the heading coupled with digested fragments (random substrings) stored in a Vectors with string elements.
  // Each sample should have equal amounts of sequencing depth. That SHOULD correspond to the number of times looped over a specific genome(?)
  // -- This means that different compositions of microbiome - i.e 1 E. coli vs 2 S. Typhimurim should be specified in the fasta file used for input.
  // The chromosome should also be circular.
  // Each fragment is at max 550 bp and lowest is 150 bp. This is nice as the tirp file will have no overlap "overflow" between the R1 and R2 columns since the smallest case will be R2 = reversed(R1). 
  // NOTE: Here, we are starting from after adapter trimming giving us our reads that are ONLY from the sample dna sequence.


  for (idx_in_fasta, (contig_name, contig_hashmap)) in &mut <Vec<(std::string::String, HashMap<&str, std::string::String>)> as Clone>::clone(&bact_entry.fasta_as_vec_hashmap).into_iter().enumerate() {
    let seq_len = contig_hashmap["contig_seq_str"].len();

    // Decide the Copy Number via poission and some valid meta knowledge
    let copy_number = if contig_hashmap["chr_name"] == "Chromosome" {1 as usize} else {poi.sample(&mut rand::rng()) as usize};
    let mut chr_name = contig_hashmap.get("chr_name").unwrap().clone();
    chr_name.push_str(&format!("{}", idx_in_fasta));

    // Float serves as "fragment per base pair". ############# CHANGE IF NEEDED OR TUNED #############
    let num_of_reads = copy_number*(0.002*(seq_len as f64)).floor() as usize;

    // Progress bar setup
    let bar = ProgressBar::new_spinner();
    bar.enable_steady_tick(Duration::from_millis(100));

    for _ in 1..=num_of_reads {  
      let start_seq_idx: usize = rand::rng().random_range(0..=seq_len);
      // Loop to ensure fragment is within length limits
      let rng_reuse = &mut rand::rng();
      let fragment_len: usize = loop {
        // Mean 400 by PI choice, stdev 50 which give 99,7 of reads being of correct length (maybe bad?)
        let fragment_len_loop: f32 = frag_len_distr.sample(rng_reuse);
        let fragment_len_loop: usize = fragment_len_loop.floor() as usize;

        // Ensuring fragment length is within bounds
        if 250 < fragment_len_loop && fragment_len_loop < 550 {
          break fragment_len_loop;
        }
      };
      // Ensure circular chromosome property for sequence
      if &start_seq_idx + &fragment_len > seq_len {
        let mut fragment = String::new();
        fragment.push_str(&contig_hashmap["contig_seq_str"][start_seq_idx..]);
        fragment.push_str(&contig_hashmap["contig_seq_str"][..(start_seq_idx+fragment_len-seq_len)]);
        format_and_write_to_tirp_line((&contig_name, &fragment), output_file, &phred_score, &cellid_hashnum);
      } else {
        let fragment = contig_hashmap["contig_seq_str"][start_seq_idx..(&start_seq_idx+&fragment_len)].to_string();
        format_and_write_to_tirp_line((&contig_name, &fragment), output_file, &phred_score, &cellid_hashnum,);
      }

  }
  bar.finish_and_clear();
  println!("Contig reads for cell#{:06} generated!", cellid_hashnum);
  

  // Write info out to the metafile
  metafile_line_writer(out_metafile, &chr_name, &cellid_hashnum, &copy_number, &num_of_reads, &bact_entry.strain_name,);
}
write_chrom_to_tirp_line(&bact_entry, output_file, &cellid_hashnum, bytes_chrom_r1);
                   }

// Writer functions
fn format_and_write_to_tirp_line(tup_contigname_fragment: (&String, &String),
                                output_file: &mut BufWriter<File>,
                                phred_score: &String,
                                cellid_hashnum: &usize,) -> () {
  // Purpose is to use the simulated fragments, reformat them to reads and write to a .tirp file.
  // --Thinking that maybe using the identifier as the cell-id?
  // Init an output file
  // Add to main()
  write!(output_file, "cell#{:06}", cellid_hashnum).expect("Problem with writing tirp data content");
    // Colidx 1 to 3 is nonsense added in written string UPDATE: I see it is 1 to 2 now.
  output_file.write_all(b"\t1").expect("Problem with writing tirp data content"); 
  output_file.write_all(b"\t1\t").expect("Problem with writing tirp data content");  
    // Make r1 and r2, cols of idx 4 and 5
  output_file.write_all(tup_contigname_fragment.1[..150].as_bytes()).expect("Problem with writing tirp data content");
  output_file.write_all(b"\t").expect("Problem with writing tirp data content"); 
  output_file.write_all(String::from_utf8(tup_contigname_fragment.1[tup_contigname_fragment.1.len()-150..tup_contigname_fragment.1.len()].chars().rev().map(|b| BASE_COMP_LUT[b as usize]).collect::<Vec<u8>>()).unwrap().as_bytes()).expect("Problem with writing tirp data content"); 
  output_file.write_all(b"\t").expect("Problem with writing tirp data content");
  // q1 and q2 (colidx 6 to 7) will point to phred_score and last col is a blankspace
  output_file.write_all(phred_score.as_bytes()).expect("Problem with writing tirp data content");
  output_file.write_all(b"\t").expect("Problem with writing tirp data content");
  output_file.write_all(phred_score.as_bytes()).expect("Problem with writing tirp data content"); 
  output_file.write_all(b"\t \n").expect("Problem with writing tirp data content"); 
  
    
}
fn write_chrom_to_tirp_line(bact_entry: &Bactdatafromfasta,
                                output_file: &mut BufWriter<File>,
                                cellid_hashnum: &usize,
                                bytes_chrom_r1: usize) -> () {
  // Purpose is to use the simulated fragments, reformat them to reads and write to a .tirp file.
  // --Thinking that maybe using the identifier as the cell-id?
  if 
  let full_seq = &bact_entry.fasta_as_vec_hashmap.clone().into_iter().filter(|(_v,hm)| hm["chr_name"] == "Chromosome").next().unwrap().1.get("contig_seq_str").unwrap().to_string();
  for chunk_idx in 0..=((full_seq.len() / bytes_chrom_r1) + 1) {
    // Quickout in case of overflow
    if chunk_idx*bytes_chrom_r1-(chunk_idx*31) >= full_seq.len(){
      break;
    }
    let seq = 
      if (chunk_idx+1)*bytes_chrom_r1-(31*chunk_idx) >= full_seq.len() {&full_seq[chunk_idx*bytes_chrom_r1-(chunk_idx*31)..]}
      else {&full_seq[chunk_idx*bytes_chrom_r1-(31*chunk_idx)..(chunk_idx+1)*bytes_chrom_r1-(chunk_idx*31)]};

    write!(output_file, "cell#{:06}_chromref", cellid_hashnum).expect("Problem with writing tirp data content");
    output_file.write_all(b"\t1").expect("Problem with writing tirp data content"); 
    output_file.write_all(b"\t1\t").expect("Problem with writing tirp data content");  
    output_file.write_all(seq.as_bytes()).expect("Problem with writing tirp data content");
    output_file.write_all(b"\t").expect("Problem with writing tirp data content");
    output_file.write_all(String::from_utf8(seq.chars().rev().map(|b| BASE_COMP_LUT[b as usize]).collect::<Vec<u8>>()).unwrap().as_bytes()).expect("Problem with writing tirp data content"); 
    output_file.write_all(b"\t").expect("Problem with writing tirp data content");
    output_file.write_all((0..seq.len()).map(|_| "F").collect::<String>().as_bytes()).expect("Problem with writing tirp data content");
    output_file.write_all(b"\t").expect("Problem with writing tirp data content");
    output_file.write_all((0..seq.len()).map(|_| "F").collect::<String>().as_bytes()).expect("Problem with writing tirp data content"); 
    output_file.write_all(b"\t \n").expect("Problem with writing tirp data content"); 
    };

}
   
fn metafile_line_writer(out_metafile: &mut BufWriter<File>,
                        chr_name: &String,
                        cellid_hashnum: &usize,
                        copy_number: &usize,
                        num_of_reads: &usize,
                        strain_name: &String,
                    ) -> () {
    let out_str_line = format!("cell#{:06}\t{}\t{}\t{}\t{}\n",cellid_hashnum, strain_name, copy_number, chr_name, num_of_reads,);
    let _ = out_metafile.write_all(out_str_line.as_bytes()).expect("Problem with writing metadata content");
}

// Source files setup functions
fn instantiate_bact(entry:DirEntry) -> Vec<(Contigfromfasta<'a>)> {
      
      let entry_path = entry.path();

      let seq_report = find_file_in_asmbly(&mut fs::read_dir(&entry_path).unwrap(), "sequence_report.jsonl");

      let asmbly_report = json::parse(&find_file_in_asmbly(&mut fs::read_dir(&entry_path).unwrap(), "assembly_data_report.jsonl")).expect("assembly json could not be parsed");

      let bactfromfasta = Bactdatafromfasta{
        fasta_as_vec_hashmap: parse_fasta(find_file_in_asmbly(&mut fs::read_dir(&entry_path).unwrap(), ".fna"), seq_report),
        strain_name: asmbly_report["checkmInfo"]["checkmMarkerSet"].to_string(),
        assembly_name: asmbly_report["currentAccession"].to_string(),
      };
bactfromfasta
}
fn find_file_in_asmbly(assembly_dir: &mut ReadDir, f_name_pat: &str) -> String {
    let res = fs::read_to_string(assembly_dir.find(|f| f.as_ref().unwrap().file_name().into_string().unwrap().ends_with(f_name_pat)).unwrap().expect("File could not be opened").path()).expect("file with ending {f_name_pat} could not be read!");
    return res;
}
fn bact_from_fasta<'a>(entry:DirEntry) -> Vec<(Contigfromfasta<'a>)> {
  // Taken from Reddit (https://www.reddit.com/r/rust/comments/r5je0y/help_parsing_a_fasta_file/) lol
    let contig_vec = fs::read_to_string(entry.path())
                              .expect("Could not read fasta file!")
                              .split(">")
                              .skip(1) // Ignores first empty split - so the empty "" that happens before the first ">"
                              .map(|s| {
                                    let contig_struct = Contigfromfasta {
                                      name: s.split("\n").next().unwrap().to_string(),
                                      seq_str: s.split("\n").skip(1).collect::<String>(),
                                      chr_name: String::new(), // Placeholder, will be filled below
                                      strain_name: String::new(), // Placeholder
                                      assembly_name: String::new(), // Placeholder
                                    };
                                    contig_struct
                                    }).collect::<Vec<(Contigfromfasta<'a>)>>();

    file_vec
}