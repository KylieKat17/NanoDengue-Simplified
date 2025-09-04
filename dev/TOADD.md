### Adding for Simplicity

- folder creation based on
- add a conf for the path where 
- just set outdir based on indir and file name
- auto set path of ref file

- update README 
  - so that there are instructions for where to put the bio-configs folder (and file)
    - should be in $HOME/bio-configs/nanoDengue.conf
  - instructions on how to update the .bashrc / .zshrc so that everything is on the PATH
  - how to run a project now



**Where to put nanoDengue.conf**

```bash
mkdir -p $HOME/bio-configs          # if it doesn’t exist
nano $HOME/bio-configs/nanoDengue.conf   # or use your favourite editor
```

Then run the pipeline (assuming `bio-toolbox/bin` is on `$PATH`):

```bash
export BIO_CONF_DIR=$HOME/bio-configs   # put this in .bashrc for permanence
nanoDengue.sh                           # auto-detects everything else
```

That’s it—feel free to add or remove lines as your environment changes!
