use predicates::prelude::*;
use snapbox::cmd::{self, Command};
use std::path::Path;

pub mod common;
use crate::common::*;

use sketchlib::sketch::multisketch::MultiSketch;

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_delete_sketches() {
        let sandbox = TestSetup::setup();

        let rfile1 = sandbox.create_rfile(&vec![
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
            "R6.fa.gz",
            "TIGR4.fa.gz",
        ]);
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .args(&["-f", rfile1])
            .arg("-v")
            .args(&["-o", "full_db"])
            .assert()
            .success();

        let rfile2 = sandbox.create_rfile(&vec![
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
        ]);
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .args(&["-f", rfile2])
            .arg("-v")
            .args(&["-o", "result_db"])
            .assert()
            .success();

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("delete")
            .arg("full_db")
            // .arg(sandbox.file_string("R6.fa.gz", TestDir::Input))
            // .arg(sandbox.file_string("TIGR4.fa.gz", TestDir::Input))
            .arg(sandbox.file_string("delete_test.txt", TestDir::Input))
            .arg("deleted_db")
            .assert()
            .success();

        let predicate_file = predicate::path::eq_file(Path::new(
            &sandbox.file_string("deleted_db.skd", TestDir::Output),
        ));
        assert_eq!(
            true,
            predicate_file.eval(Path::new(
                &sandbox.file_string("result_db.skd", TestDir::Output)
            )),
            "Merged sketch data does not match"
        );

        // Check .skm the same
        let merged_sketch: MultiSketch =
            MultiSketch::load_metadata(&sandbox.file_string("deleted_db", TestDir::Output))
                .expect("Failed to load output merged sketch");
        let expected_sketch =
            MultiSketch::load_metadata(&sandbox.file_string("result_db", TestDir::Output))
                .expect("Failed to load expected merged sketch");
        assert_eq!(
            merged_sketch, expected_sketch,
            "Merged sketch metadata does not match"
        );
    }

    #[test]
    /// Regression test: deleting a *non-trailing* sample (here, the second of
    /// four) used to corrupt the output `.skd`. `remove_genomes` was called
    /// after `remove_metadata` had already truncated `sketch_metadata`, so
    /// `indices_to_keep` was built from the wrong (post-truncation) length
    /// and silently dropped trailing on-disk records that should have been
    /// kept. The bug was masked by the original `test_delete_sketches` test,
    /// which only ever deletes a trailing block.
    fn test_delete_non_trailing_sample() {
        let sandbox = TestSetup::setup();

        let rfile1 = sandbox.create_rfile(&vec![
            "14412_3#82.contigs_velvet.fa.gz",
            "14412_3#84.contigs_velvet.fa.gz",
            "R6.fa.gz",
            "TIGR4.fa.gz",
        ]);
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .args(&["-f", rfile1])
            .arg("-v")
            .args(&["-o", "full_db"])
            .assert()
            .success();

        let rfile2 = sandbox.create_rfile(&vec![
            "14412_3#82.contigs_velvet.fa.gz",
            "R6.fa.gz",
            "TIGR4.fa.gz",
        ]);
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .args(&["-f", rfile2])
            .arg("-v")
            .args(&["-o", "result_db"])
            .assert()
            .success();

        std::fs::write(
            sandbox.file_string("delete_middle.txt", TestDir::Output),
            "14412_3#84.contigs_velvet.fa.gz\n",
        )
        .expect("Failed to write delete list");

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("delete")
            .arg("full_db")
            .arg("delete_middle.txt")
            .arg("deleted_db")
            .assert()
            .success();

        let predicate_file = predicate::path::eq_file(Path::new(
            &sandbox.file_string("deleted_db.skd", TestDir::Output),
        ));
        assert!(
            predicate_file.eval(Path::new(
                &sandbox.file_string("result_db.skd", TestDir::Output)
            )),
            "Deleting a non-trailing sample produced a .skd that doesn't match a fresh sketch of the remaining samples"
        );

        let deleted_sketch: MultiSketch =
            MultiSketch::load_metadata(&sandbox.file_string("deleted_db", TestDir::Output))
                .expect("Failed to load deleted sketch");
        assert_eq!(deleted_sketch.number_samples_loaded(), 3);
        assert_eq!(
            deleted_sketch.get_sample_index("R6.fa.gz"),
            Some(1),
            "name_map should map retained samples to their new compacted position"
        );
        assert_eq!(deleted_sketch.get_sample_index("TIGR4.fa.gz"), Some(2));
        assert_eq!(
            deleted_sketch.get_sample_index("14412_3#84.contigs_velvet.fa.gz"),
            None,
            "name_map should not retain the deleted sample"
        );
    }

    #[test]
    /// Regression test: `remove_metadata` never used to update `name_map`, so
    /// a deleted sample's name stayed in it. Merging two databases each
    /// produced by deleting down to a disjoint single sample would then
    /// spuriously panic with "<name> appears in both databases", even though
    /// that sample was never actually present on the other side.
    fn test_merge_after_delete_does_not_see_stale_names() {
        let sandbox = TestSetup::setup();

        let rfile = sandbox.create_rfile(&vec!["R6.fa.gz", "TIGR4.fa.gz"]);
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("sketch")
            .args(&["--k-vals", "17"])
            .args(&["-f", rfile])
            .arg("-v")
            .args(&["-o", "pair_db"])
            .assert()
            .success();

        std::fs::write(
            sandbox.file_string("remove_tigr4.txt", TestDir::Output),
            "TIGR4.fa.gz\n",
        )
        .expect("Failed to write delete list");
        std::fs::write(
            sandbox.file_string("remove_r6.txt", TestDir::Output),
            "R6.fa.gz\n",
        )
        .expect("Failed to write delete list");

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("delete")
            .arg("pair_db")
            .arg("remove_tigr4.txt")
            .arg("r6_only")
            .assert()
            .success();
        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("delete")
            .arg("pair_db")
            .arg("remove_r6.txt")
            .arg("tigr4_only")
            .assert()
            .success();

        Command::new(cmd::cargo_bin!("sketchlib"))
            .current_dir(sandbox.get_wd())
            .arg("merge")
            .arg("r6_only")
            .arg("tigr4_only")
            .args(&["-o", "remerged"])
            .assert()
            .success();

        let merged: MultiSketch =
            MultiSketch::load_metadata(&sandbox.file_string("remerged", TestDir::Output))
                .expect("Failed to load remerged sketch");
        assert_eq!(merged.number_samples_loaded(), 2);
        assert_eq!(merged.get_sample_index("R6.fa.gz"), Some(0));
        assert_eq!(merged.get_sample_index("TIGR4.fa.gz"), Some(1));
    }
}
