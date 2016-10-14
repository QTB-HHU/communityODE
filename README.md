## ODE model for P. tricornutum associated microbial community

Please refer to the project wiki for more info

### Installation

This code requires some external open-source packages.

We use virtual environment, if you don't have it please install like one of:

```bash
apt-get install python-virtualenv
dnf install python-virtualenv
pip install virtualenv
```

Generate a virtual environment for python2 or python3 like:

```bash
virtualenv -p /yourlocalpathto/pythonX yourprojectfolder
```

e.g.

```bash
virtualenv -p /usr/bin/python2.7 $HOME/diatom-comm/code/python
```


Then activate virtualenv and install the requirements:

```bash
cd yourprojectfolder
source bin/activate
pip install -r requirements.txt
```

You are ready to run! When done, leave the virtualenv like:

```bash
deactivate
```

N.B. when freshly creating virtualenv on fedora, installing via requirements all at once does not work... this worked:
pip install numpy
pip install numexpr bottleneck
pip install six cycler
pip install pytz
pip install -r requirements.txt
