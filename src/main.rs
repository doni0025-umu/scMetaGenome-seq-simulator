//Program should read a genome file, substring it into reads and finally write it all to a "tirp" (tab-indexed-read-paired) file.

use std::collections::HashMap;
use std::fs;
use std::env;
use std::fs::File;
use std::fs::ReadDir;
use rand::{rng, Rng};
use rand_distr::{Normal, Distribution, Poisson};
use std::io::Write;
use json;
use std::fs::DirEntry;

struct Bactdatafromfasta <'a>{
  fasta_as_vec_hashmap: Vec<(String,HashMap<&'a str,String>)>,
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
    let poi = Poisson::new(13.0).unwrap();

    // Preamble for write-to-tirp
    let out_path = &args[1];
    let mut output_file  = std::fs::OpenOptions::new().write(true).truncate(true).create(true).open(out_path).expect("This path could NOT be used as output path. Maybe dir does not exist?");
    let phred_score_prea = (0..150).map(|_| "F").collect::<String>();
    let base_comp_prea = HashMap::from([('A', 'T'), ('T','A'),('G', 'C'), ('C','G'), ('N', 'N'),]);

    // Preamble for metafile_line_writer
    let out_metafile_path = &args[2];
    let out_metafile = std::fs::OpenOptions::new().write(true).truncate(true).create(true).open(out_metafile_path).expect("This path could NOT be used as Metafile output path. Maybe dir does not exist?");

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
        read_simulator(&bact_entry, min_len, max_len, fragm_per_bp, frag_len_distr, &mut output_file, &base_comp_prea, &phred_score_prea, cellid_hashnum.id_counter, &out_metafile, poi);
        cellid_hashnum.count()
      }
    }
    Ok(())
    }


// Simulation function
fn read_simulator(bact_entry: &Bactdatafromfasta,
                    min_len: usize,
                    max_len: usize, 
                    frag_per_bp: f64, 
                    frag_len_distr: Normal<f32>,
                    output_file: &mut File, 
                    base_comp: &HashMap<char,char>,
                    phred_score: &String,
                    cellid_hashnum: usize,
                    out_metafile: &File,
                    poi: Poisson<f64>,
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

    let num_of_reads = copy_number*(frag_per_bp*(seq_len as f64)).floor() as usize;
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
        fragment.push_str(&contig_hashmap["contig_seq_str"][start_seq_idx..]);
        fragment.push_str(&contig_hashmap["contig_seq_str"][..(start_seq_idx+fragment_len-seq_len)]);
        format_and_write_to_tirp_line((&contig_name, &fragment), &output_file, &base_comp, &phred_score, &cellid_hashnum);
      } else {
        let fragment = contig_hashmap["contig_seq_str"][start_seq_idx..(&start_seq_idx+&fragment_len)].to_string();
        format_and_write_to_tirp_line((&contig_name, &fragment), &output_file, &base_comp, &phred_score, &cellid_hashnum,);
      }

  }
  // Write info out to the metafile
  metafile_line_writer(&out_metafile, &chr_name, &cellid_hashnum, &copy_number, &num_of_reads, &bact_entry.strain_name,);
}
write_chrom_to_tirp_line(&bact_entry, output_file, &cellid_hashnum,)
                   }

// Writer functions
fn format_and_write_to_tirp_line(tup_contigname_fragment: (&String, &String),
                                mut output_file: &File, 
                                base_comp: &HashMap<char,char>, 
                                phred_score: &String,
                                cellid_hashnum: &usize,) -> () {
  // Purpose is to use the simulated fragments, reformat them to reads and write to a .tirp file.
  // --Thinking that maybe using the identifier as the cell-id?
  // Init an output file
  // Add to main()
  let fragment = tup_contigname_fragment.1;
  // Colidx 1 to 3 is nonsense added in written string UPDATE: I see it is 1 to 2 now.
  // Make r1 and r2, cols of idx 4 and 5
  let r1 = &fragment[..150];
  // Reverses fragment and complements the 150 bases via a HashMap
  let r2 = &fragment[fragment.len()-150..fragment.len()].chars().rev().map(|b|{base_comp.get(&b).copied().unwrap_or('N')}).collect::<String>();

  // q1 and q2 (colidx 6 to 7) will point to phred_score and last col is a blankspace
  let out_str_line = format!("cell#{:06}\t1\t1\t{}\t{}\t{}\t{}\t \n", cellid_hashnum, r1, r2, &phred_score, &phred_score);
  let _ = output_file.write_all(out_str_line.as_bytes()).expect("Problem with writing tirp data content");
    
}
fn write_chrom_to_tirp_line(bact_entry: &Bactdatafromfasta,
                                output_file: &mut File,
                                cellid_hashnum: &usize,) -> () {
  // Purpose is to use the simulated fragments, reformat them to reads and write to a .tirp file.
  // --Thinking that maybe using the identifier as the cell-id?
  let seq = &bact_entry.fasta_as_vec_hashmap.clone().into_iter().filter(|(_v,hm)| hm["chr_name"] == "Chromosome").next().unwrap().1.get("contig_seq_str").unwrap().to_string();
  let out_str_line = format!("cell#{:06}_chromref\t1\t1\t{}\t\t{}\t\t \n", cellid_hashnum, seq, (0..seq.len()).map(|_| "F").collect::<String>());
  let _ = output_file.write_all(out_str_line.as_bytes()).expect("Problem with writing tirp data content");
          
        }
fn metafile_line_writer(mut out_metafile: &File,
                        chr_name: &String,
                        cellid_hashnum: &usize,
                        copy_number: &usize,
                        num_of_reads: &usize,
                        strain_name: &String,
                    ) -> () {
/*  cellID	strain	copyNumber	contigName	readCount
    0001	salmonella	1	main	    12
    0001	salmonella	2	plasmid1	1234
    0001	salmonella	3	plasmid2	234 */
    let out_str_line = format!("cell#{:06}\t{}\t{}\t{}\t{}\n",cellid_hashnum, strain_name, copy_number, chr_name, num_of_reads,);
    let _ = out_metafile.write_all(out_str_line.as_bytes()).expect("Problem with writing metadata content");
    println!("Metaline for cell#{:06} written!", cellid_hashnum);
}

// Source files setup functions
fn instantiate_bact(entry:DirEntry) -> Bactdatafromfasta<'static> {
      
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
fn parse_fasta<'a>(file: String, seq_report: String) -> Vec<(String,HashMap<&'a str,String>)> {
  // Taken from Reddit (https://www.reddit.com/r/rust/comments/r5je0y/help_parsing_a_fasta_file/) lol
    let file_vec = file.split(">")         
      .skip(1) // Ignores first empty split - so the empty "" that happens before the first ">"
      .map(|s| {
        let mut lines = s.split("\n");
        let contig_name = lines.next().unwrap().to_string();
        let mut contig_hashmap = HashMap::from([("contig_seq_str",lines.collect::<String>()), ("refSeqAccession", contig_name.split_whitespace().next().unwrap().to_string())]);
        let json_data = json::parse(&seq_report.lines().filter(|line| line.contains(contig_hashmap.get("refSeqAccession").unwrap())).next().unwrap()).expect(&format!("Not able to parse seqreport for contig ****{}****.", contig_name));
        contig_hashmap.insert("chr_name", json_data["assignedMoleculeLocationType"].to_string(),);

        (contig_name, contig_hashmap)
      })
      .collect::<Vec<(String,HashMap<&str,String>)>>();
    file_vec
}