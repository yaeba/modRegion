# modRegion
 A tool to quickly access, query and merge per-read base modification results from different software tools

## Example Usage
```bash
time ./modRegion.R extract -r "NC_001144.5:430000-470000" -r "NC_001144.5:400000-420000" \
  -o test/small_deepsignal.extracted.tsv.gz \
  --deepsignal test/small_deepsignal.tsv.gz \
  -p test/small_deepsignal.plots.pdf
```

```bash
time ./modRegion.R extract -r "OCVW01000001.1:40,000-80,000"  \
  -o test/small_nanotom.extracted.tsv.gz \
  --nanopolish nanopolish:test/small_nanopolish.tsv.gz \
  --tombo tombo:test/small_tombo.tsv.gz \
  -p test/small_nanotom.plots.pdf
```
