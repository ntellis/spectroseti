spectroseti
======
Tools for searching for narrowband optical emission in spectroscopic data

About the code
--------------------------
``spectroseti`` contains all tools you will need to search spectroscopic data
for laser lines. This package was developed by Nate Tellis for the Breakthrough
Listen SETI program.

WIP



Installation
--------------------------

You can easily configure a local instance of python with all of the package dependencies, but if you do not want to do this, you can use ``miniconda``.

Currently ``spectroseti`` supports only python 2.7.x


##### To install with ``miniconda``

```commandline
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod +x Miniconda2-latest-Linux-x86_64.sh
./Miniconda2-latest-Linux-x86_64.sh
```

This should have added miniconda to your PATH, but if it failed, add this line to your .bashrc:
```commandline
export PATH="/INSTALLDIR/miniconda2/bin:$PATH"
```
where `INSTALLDIR` is the directory to which miniconda was installed (by default /home/your_user_name/)



Dependencies:
numpy
matplotlib
pandas
pathos (and dill)
tqdm
scipy
astropy
random
sklearn


WIP



Acknowledgments


This package makes use of the Pathos library for expanded multiprocessing in Python

M.M. McKerns, L. Strand, T. Sullivan, A. Fang, M.A.G. Aivazis,
"Building a framework for predictive science", Proceedings of
the 10th Python in Science Conference, 2011;
http://arxiv.org/pdf/1202.1056

Michael McKerns and Michael Aivazis,
"pathos: a framework for heterogeneous computing", 2010- ;
http://trac.mystic.cacr.caltech.edu/project/pathos
