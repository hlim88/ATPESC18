## Setup on a cooley login node

```
soft add +anaconda3-4.0.0

# generate config file
# ~/.jupyter/jupyter_notebook_config.py
jupyter notebook --generate-config

# Generate an encrypted password 
# Start a python shell and run
from notebook.auth import passwd; passwd()

enter and confirm password and you
should get something like
'sha1:abdnadnabadadbnandabdandnabdadadadadada`
```

Now edit
```
~/.jupyter/jupyter_notebook_config.py
```

and enter the encrypted password in line 201

```
c.NotebookApp.password = 'sha1:abdnadnabadadbnandabdandnabdadadadadada`

# Enable only local network connections
# line 155
c.NotebookApp.ip = 'localhost'

# no browser
# line 192
c.NotebookApp.open_browser = False

# line 204 - select a port in the range 8000 - 8008
# The port the notebook server will listen on.
c.NotebookApp.port = 8000
```

---
**^ ^ ^ NOTE THAT YOU MUST UNCOMMENT THE LINES
ON THE FILE ABOVE FOR IT TO WORK ^ ^ ^**

---

## Interactive session on cooley

From login node launch an interactive session on a single node. Submit job to *training* queue and use ATPESC2018 allocation

```
qsub -n 1 -I -t 60 -q training -A ATPESC2018
```

When you get the interactive prompt make a note of your node name (i.e cc008)

```
# set environment
soft add +anaconda3-4.0.0

# Launch jupyter notebook
jupyter-notebook

# This will launch jupyter server and start listening on
http://cc008.cooley.pub.alcf.anl.gov:8000/
```

## From your laptop or desktop

```
# From another terminal create an ssh tunnel to your
# compute node passing through login node
# This is to make sure your traffic (and password)
# are encrypted

# Note we use cc008 as an example here, but you must replace it
# with the name of the node you got in the step above.
# We are also using port 9000 as an example here and you
# are free to change it
ssh -l your_username -L 9000:localhost:9000 cooley.alcf.anl.gov  ssh -L 9000:localhost:8000 -N cc008

# point browser to your local end of the tunnel
http://localhost:9000/
```
