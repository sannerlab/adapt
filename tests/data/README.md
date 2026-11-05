# Entries from the RCSB used by ADapt's tests

Downloaded from `https://files.rcsb.org/download/<ID>.pdb` and gzipped.

| File | Downloaded | md5 of the .gz | Used for |
|---|---|---|---|
| 4ins.pdb.gz | 2026-11-05 | eda05fb061903f1d48826061d56cc304 | insulin hexamer (Biomol 3): zinc ions and waters on the 3-fold axis are written once, not once per copy (test_copies) |
| 1plc.pdb.gz | 2026-11-05 | f4c0e787b663938de4d7d7f157c09013 | plastocyanin copper replaced by zinc from the tree, and undone (test_choices) |
| 1ddy.pdb.gz | 2026-11-05 | 77fe6fb56afa7058ab40c0dc8337ac60 | B12 (cobalt) left out of a receptor, reported in the save record (test_choices) |
