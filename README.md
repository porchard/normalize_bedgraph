# Normalize bedgraph
Given a bedgraph file, normalize it to a given signal/base or to a given
number of reads.

## Usage:
Pass the `-h` flag for a list of arguments. Two common use cases would be:
```bash
python normalize_bedgraph.py --to-mean-signal 1.0 unnormalized.bdg > normalized.bdg
python normalize_bedgraph.py --to-number-reads 10000000 unnormalized.bdg > normalized.bdg
```
