.. _installation:

Installation 
------------

For now, you should check out the anndata branch of Monod directly from Github, as in the demo notebook.

 .. code-block:: bash

  git clone https://github.com/pachterlab/monod.git
  git checkout anndata

Install any dependencies:

 .. code-block:: bash

  pip install -q ./
  pip install -q -r /content/monod/src/monod.egg-info/requires.txt
  pip install scanpy

TODO: make anndata version pip installable.

..
 (NB not valid currently - we need to update the Monod version. For now ignore this section and use git pull instead as in demo notebook)
 
 To use *Monod*, install it from `pip` (TODO update version):
 
 .. code-block:: console
 
  pip install monod
  
 To use it in your code, import the package components:
 
 .. code-block:: python
 
  import monod
  from monod import *
  from monod.extract_data import *
  from monod.cme_toolbox import CMEModel
 from monod import inference, mminference
 from monod.analysis import *
