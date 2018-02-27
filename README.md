# EMEP/MSC-W Model Unofficial User’s Guide
[![Documentation Status](https://readthedocs.org/projects/emep-ctm/badge/?version=latest)](http://emep-ctm.readthedocs.io/en/latest/?badge=latest)

This is part of the [Open Source EMEP/MSC-W model][emep-ctm] documentation.

[emep-ctm]: https://github.com/metno/emep-ctm

## build environment
The [user guide][] is build with with [sphinx][],
installed in a [conda environment][conda].

[user guide]: User_Guide.pdf
[sphinx]: http://www.sphinx-doc.org
[conda]: http://conda.pydata.org

### LaTeX package
```bash
sudo apt-get install latexmk
```

### conda setup
```bash
# create env
conda create --name sphinx sphinx

# activate env
source activate sphinx

# add modules
# conda install sphinx-autobuild
pip install sphinxcontrib-bibtex
```

### run sphinx
```bash
# update build/latex/emep-ctm.pdf
make latexpdf   

# hard link pdf
ln -s build/latex/emep-ctm.pdf User_Guide.pdf
```
