# ODE model for P. tricornutum associated microbial community

Pre-print doi: [http://dx.doi.org/10.1101/077768](http://dx.doi.org/10.1101/077768)


## Content

This model defines a Ordinary Differential Equations system to describe the dynamic
evolution of an ecosystem of five organisms (one diatom and four bacteria). 
[//]: # (Also provided is a genetic algorithm code used to define the model parameters.)

The code is made available under the GNU General Public License (see LICENSE) at no warranty.

Macros to run are available in:

```bash
macros/examples/
macros/diatom-comm/
```

[//]: # (Code documentation is available in:)


## Installation

This code requires some external open-source packages.

We use virtual environment, if you don't have it please install like one of:

```bash
apt-get install python-virtualenv
dnf install python-virtualenv
pip install virtualenv
```

Generate a virtual environment for python2 or python3 like:

```bash
virtualenv -p /yourlocalpathto/pythonX yourprojectfolder/venv
```

e.g.

```bash
virtualenv -p /usr/bin/python2.7 $HOME/diatom-comm/venv
```


Then activate virtualenv and install the requirements:

```bash
cd yourprojectfolder
source venv/bin/activate
pip install -r code/python/requirements.txt
```

You are ready to run! When done, leave the virtualenv like:

```bash
deactivate
```
