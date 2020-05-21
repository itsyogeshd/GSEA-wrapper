# Wrapper script for running GSEA via command line

```bash
python gseawrap.py --help
usage: Wrapper for GSEA. [-h] [--gct GCT] [--cls CLS] [--gmt GMT]
                         [--chip CHIP] [--projectname PROJECTNAME]
                         [--nplots NPLOTS] [--nperms NPERMS]
                         [--metric {Signal2Noise,tTest,log2_Ratio_of_Classes}]
                         [--ispreranked] [--rnk [RNK [RNK ...]]]
                         [--log-level {NOTSET,DEBUG,INFO,WARNING,CRITICAL,ERROR,CRITICAL}]

optional arguments:
  -h, --help            show this help message and exit
  --gct GCT             Gene expression matrix in .gct format.
  --cls CLS             Path to class file in .cls format.
  --gmt GMT             Path to geneset file in .gmt format.
  --chip CHIP           Path to annotation file in .chip format.
  --projectname PROJECTNAME
                        Prefix for project with absolute path.
  --nplots NPLOTS       Number of plots to generate.
  --nperms NPERMS       Number of permutations.
  --metric {Signal2Noise,tTest,log2_Ratio_of_Classes}
                        Gene ranking method to choose (If number of samples in
                        each group is more than 3, choose "Signal2Noise") (If
                        number of samples in each group is less than 3, choose
                        either "tTest" or "log2_Ratio_of_Classes")
  --ispreranked         Analysis is preranked GSEA.
  --rnk [RNK [RNK ...]]
                        Preranked genelist in .rnk format.
  --log-level {NOTSET,DEBUG,INFO,WARNING,CRITICAL,ERROR,CRITICAL}
                        Log level
```
