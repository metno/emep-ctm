# EMEP/MSC-W Model Unofficial User’s Guide

This is part of the [Open Source EMEP/MSC-W model][emep-ctm] documentation.

[emep-ctm]: https://github.com/metno/emep-ctm

## build environment
The [user guide][] is build with with [sphinx][],
installed in a [conda environment][conda].

[user guide]: User_Guide.pdf
[sphinx]: http://www.sphinx-doc.org
[conda]: http://conda.pydata.org

### conda setup
```bash
# create env
conda create --name sphinx sphinx

# activate env
source activate sphinx

# add modules
# conda install sphinx-autobuild
pip install sphinxcontrib-bibtexls
```

### run sphinx
```bash
# update build/latex/emep-ctm.pdf
make latexpdf   

# hard link pdf
ln -s build/latex/emep-ctm.pdf User_Guide.pdf
```
