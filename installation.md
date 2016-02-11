---
layout: page
title: Installation  
subtitle: Installing RapClust 
---

## Dependencies
----------------

RapClust depends on the following external programs (to be available in the environment where it runs):

  1. The [MCL](http://micans.org/mcl/) clustering tool
  2. The [Sailfish](https://github.com/kingsfordgroup/sailfish) (or [Salmon](https://github.com/COMBINE-lab/salmon)) quantification tool<sup id="a1">[1](#f1)</sup>

Further, it depends on the following Python packages:
  
  1. [Click](http://click.pocoo.org/5/)
  2. [PyYAML](https://pypi.python.org/pypi/PyYAML)
  3. [Pandas](http://pandas.pydata.org/)
  4. [NumPy](http://www.numpy.org/)

However, you should be able to install rapclust via `pip` and have these dependencies installed automatically.  To install RapClust via pip, you can use:

## Installation via pip
-----------------------

```
> pip install rapclust
```

You should now have a `RapClust` executable in your path.  You can test this with the following command:

```
> RapClust --help
```

You should see the following output:

```
Usage: RapClust [OPTIONS]

Options:
  --config TEXT  Config file describing the experimental setup
  --help         Show this message and exit.
```

If you see this, then RapClust has been successfully installed!

