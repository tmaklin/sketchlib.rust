use snapbox::cmd::{self, Command};

pub mod common;
use crate::common::*;

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn sketch_fasta() {
        let sandbox = TestSetup::setup();

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("assembly")
            .args(["-v", "-k", "31"])
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("assembly.skm"));
        assert_eq!(true, sandbox.file_exists("assembly.skd"));

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("assembly")
            .assert()
            .stdout_eq(sandbox.snapbox_file("assembly_sketch_info.stdout", TestDir::Correct));

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("--sample-info")
            .arg("assembly")
            .arg("-v")
            .assert()
            .stdout_eq(sandbox.snapbox_file("assembly_sketch_full_info.stdout", TestDir::Correct));
    }

    #[test]
    fn sketch_fastq() {
        let sandbox = TestSetup::setup();

        // Create a fastq rfile in the tmp dir
        let rfile_name = sandbox.create_fastq_rfile("test");
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-f")
            .arg(rfile_name)
            .arg("-o")
            .arg("reads")
            .args(["--min-count", "2", "-v", "-k", "9", "--min-qual", "2"])
            .assert()
            .success();

        assert_eq!(true, sandbox.file_exists("reads.skm"));
        assert_eq!(true, sandbox.file_exists("reads.skd"));

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("reads")
            .assert()
            .stdout_eq(sandbox.snapbox_file("read_sketch_info.stdout", TestDir::Correct));

        // FASTQ read length is estimated from minima normalized over the full u64 hash space.
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("info")
            .arg("--sample-info")
            .arg("reads")
            .arg("-v")
            .assert()
            .stdout_eq(sandbox.snapbox_file("read_sketch_full_info.stdout", TestDir::Correct));
    }

    #[test]
    fn sketch_fastq_bad() {
        let sandbox = TestSetup::setup();

        let rfile_name_bad = sandbox.create_bad_fastq_rfile("test");

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-f")
            .arg(rfile_name_bad)
            .arg("-o")
            .arg("readsbad")
            .args(["--min-count", "2", "-v", "-k", "9", "--min-qual", "2"])
            .assert()
            .failure();
    }

    #[test]
    fn sketch_aas() {
        let sandbox = TestSetup::setup();

        sandbox.copy_input_file_to_wd("test_aa_sequence.fa");

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(["-o", "aastest"])
            .args(["--seq-type", "aa"])
            .args(["--min-count", "2", "-v", "--k-vals", "9", "--min-qual", "2"])
            .arg("./test_aa_sequence.fa")
            .assert()
            .success();

        // check level 2
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(["-o", "aastest"])
            .args(["--seq-type", "aa"])
            .args(["--level", "level2"])
            .args(["--min-count", "2", "-v", "--k-vals", "9", "--min-qual", "2"])
            .arg("./test_aa_sequence.fa")
            .assert()
            .success();

        // check level 3
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(["-o", "aastest"])
            .args(["--seq-type", "aa"])
            .args(["--level", "level3"])
            .args(["--min-count", "2", "-v", "--k-vals", "9", "--min-qual", "2"])
            .arg("./test_aa_sequence.fa")
            .assert()
            .success();
    }

    #[test]
    /// The MultiSketch (.skd/.skm) format changed as of v0.4, so databases
    /// created with older versions (this fixture is v0.1.3) must be refused
    /// rather than silently loaded.
    fn legacy_database_rejected() {
        use sketchlib::sketch::multisketch::MultiSketch;
        let sandbox = TestSetup::setup();

        let err = MultiSketch::load_metadata(&sandbox.file_string("legacy_db", TestDir::Input))
            .expect_err("Loading a pre-v0.4 sketch file should fail");
        assert!(
            err.to_string().contains("Incompatible sketch file version"),
            "Unexpected error loading legacy sketch: {err}"
        );
    }

    #[test]
    fn load_matches_load_metadata_then_read_sketch_data() {
        use sketchlib::sketch::multisketch::MultiSketch;
        let sandbox = TestSetup::setup();

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .arg("-o")
            .arg("loadtest")
            .args(["-k", "31"])
            .arg(sandbox.file_string("14412_3#82.contigs_velvet.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("14412_3#84.contigs_velvet.fa.gz", TestDir::Input))
            .assert()
            .success();

        let prefix = sandbox.file_string("loadtest", TestDir::Output);

        let mut separate =
            MultiSketch::load_metadata(&prefix).expect("failed to load metadata separately");
        separate.read_sketch_data(&prefix);

        let combined = MultiSketch::load(&prefix).expect("failed to load combined");

        assert_eq!(
            separate, combined,
            "MultiSketch::load should match load_metadata + read_sketch_data"
        );

        // sketch_metadata() should expose the same names as the indexed accessors
        let names: Vec<&str> = combined
            .sketch_metadata()
            .iter()
            .map(|s| s.name())
            .collect();
        assert_eq!(names.len(), combined.number_samples_loaded());
        for (idx, name) in names.iter().enumerate() {
            assert_eq!(*name, combined.sketch_name(idx));
        }
    }
}
